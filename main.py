import os
import subprocess
import urllib.request
from collections import defaultdict


def read_fasta(filepath):
    """Читает FASTA файл и возвращает список последовательностей."""
    seqs = []
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    seqs.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            seqs.append("".join(current_seq))
    return seqs


def read_fastq(filepath):
    """Читает FASTQ файл и возвращает список последовательностей."""
    seqs = []
    with open(filepath, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # строка '+'
            f.readline()  # строка качества
            seqs.append(seq)
    return seqs


class DeBruijnGraph:
    def __init__(self, k):
        self.k = k
        self.edge_cov = defaultdict(int)
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)
        self.adj = defaultdict(list)
        self.nodes = set()
        self.compressed_edges = []

    def add_sequence(self, seq):
        """Добавляет последовательность в граф, извлекая (k+1)-меры."""
        seq = seq.upper().replace('N', '')
        seq_len = len(seq)
        for i in range(seq_len - self.k):
            k1_mer = seq[i:i + self.k + 1]
            self.edge_cov[k1_mer] += 1

    def build_structure(self):
        """
        Строит списки смежности и степени вершин на основе накопленных рёбер.
        Вызывается один раз после добавления всех последовательностей или после очистки.
        """
        self.adj = defaultdict(list)
        self.in_degree = defaultdict(int)
        self.out_degree = defaultdict(int)
        self.nodes = set()

        for k1_mer in self.edge_cov:
            u = k1_mer[:-1]
            v = k1_mer[1:]
            self.nodes.add(u)
            self.nodes.add(v)
            self.adj[u].append(k1_mer)
            self.out_degree[u] += 1
            self.in_degree[v] += 1

    def compress(self):
        """Сжимает линейные пути графа. Вершины остаются k-мерами, рёбра содержат полные последовательности."""
        branching_nodes = set()
        for node in self.nodes:
            if self.in_degree[node] != 1 or self.out_degree[node] != 1:
                branching_nodes.add(node)

        # Обработка изолированных циклов
        visited = set()
        for node in self.nodes:
            if node not in visited:
                curr = node
                path = []
                while (curr not in visited and curr in self.adj and
                       len(self.adj[curr]) == 1 and self.in_degree[curr] == 1):
                    visited.add(curr)
                    path.append(curr)
                    curr = self.adj[curr][0][1:]
                if curr in path:
                    branching_nodes.add(curr)

        self.compressed_edges = []
        visited_k1_mers = set()

        for start_node in branching_nodes:
            for k1_mer in self.adj[start_node]:
                if k1_mer in visited_k1_mers:
                    continue

                curr_k1 = k1_mer
                seq = curr_k1
                cov = self.edge_cov[curr_k1]
                length = 1
                visited_k1_mers.add(curr_k1)
                curr_v = curr_k1[1:]

                while curr_v not in branching_nodes:
                    if self.out_degree[curr_v] == 0:
                        break
                    next_k1 = self.adj[curr_v][0]
                    if next_k1 in visited_k1_mers:
                        break
                    visited_k1_mers.add(next_k1)
                    seq += next_k1[-1]
                    cov += self.edge_cov[next_k1]
                    length += 1
                    curr_k1 = next_k1
                    curr_v = curr_k1[1:]

                self.compressed_edges.append({
                    'u': k1_mer[:-1],
                    'v': curr_v,
                    'seq': seq,
                    'cov': cov,
                    'len': length
                })

    def clean_uncompressed(self, min_cov=2, max_tip_len_k1=5):
        """Очищает граф: удаляет рёбра с низким покрытием и короткие тупиковые ветви."""
        # 1. Удаление рёбер с низким покрытием
        to_remove = [k1 for k1, cov in self.edge_cov.items() if cov < min_cov]
        for k1 in to_remove:
            del self.edge_cov[k1]

        # Перестроение структуры после удаления по покрытию
        self.build_structure()

        # 2. Удаление коротких тупиковых ветвей
        dead_ends = []
        for node in self.nodes:
            if (self.in_degree[node] == 0 and self.out_degree[node] == 1) or \
                    (self.out_degree[node] == 0 and self.in_degree[node] == 1):
                dead_ends.append(node)

        tips_to_remove = []

        # Словарь для обратного поиска (от стока к источнику)
        rev_adj = defaultdict(list)
        for u, k1_list in self.adj.items():
            for k1 in k1_list:
                rev_adj[k1[1:]].append((u, k1))

        for start_node in dead_ends:
            if self.in_degree[start_node] == 0 and self.out_degree[start_node] == 1:
                # Прямой тупик (источник)
                curr = start_node
                path_k1s = []
                while self.out_degree[curr] == 1 and self.in_degree[curr] <= 1:
                    k1 = self.adj[curr][0]
                    path_k1s.append(k1)
                    curr = k1[1:]
                    if self.in_degree[curr] > 1 or self.out_degree[curr] > 1:
                        break
                if len(path_k1s) <= max_tip_len_k1:
                    tips_to_remove.extend(path_k1s)

            elif self.out_degree[start_node] == 0 and self.in_degree[start_node] == 1:
                # Обратный тупик (сток)
                curr = start_node
                path_k1s = []
                while self.in_degree[curr] == 1 and self.out_degree[curr] <= 1:
                    u, k1 = rev_adj[curr][0]
                    path_k1s.append(k1)
                    curr = u
                    if self.in_degree[curr] > 1 or self.out_degree[curr] > 1:
                        break
                if len(path_k1s) <= max_tip_len_k1:
                    tips_to_remove.extend(path_k1s)

        for k1 in set(tips_to_remove):
            if k1 in self.edge_cov:
                del self.edge_cov[k1]

        # Финальное перестроение и сжатие
        self.build_structure()
        self.compress()

    def save_contigs(self, filepath):
        """Сохраняет сжатые линейные участки в формате FASTA."""
        with open(filepath, 'w') as f:
            for i, edge in enumerate(self.compressed_edges):
                f.write(f">contig_{i} len_k1={edge['len']} cov={edge['cov']}\n")
                seq = edge['seq']
                for j in range(0, len(seq), 80):
                    f.write(seq[j:j + 80] + "\n")

    def save_gfa(self, filepath):
        """Сохраняет граф в формате GFA1."""
        with open(filepath, 'w') as f:
            f.write("H\tVN:Z:1.0\n")
            nodes = set()
            for edge in self.compressed_edges:
                nodes.add(edge['u'])
                nodes.add(edge['v'])

            for node in sorted(list(nodes)):
                f.write(f"S\t{node}\t*\n")

            for edge in self.compressed_edges:
                u = edge['u']
                v = edge['v']
                seq = edge['seq']
                cov = edge['cov']
                f.write(f"L\t{u}\t+\t{v}\t+\t*\tse:Z:{seq}\tcv:i:{cov}\n")

    def get_stats(self):
        return f"Вершин (k-меров): {len(self.nodes)}, Уникальных рёбер (k+1-меров): {len(self.edge_cov)}, Сжатых рёбер: {len(self.compressed_edges)}"


def main():
    ref_file = "ecoli_1k.fna"
    reads_file = "ecoli_reads.fastq"
    bandage_exc = "./../Bandage"
    data_dir = "./data/"


    # Автоматическая загрузка файлов, если их нет
    if not os.path.exists(ref_file):
        print("Скачивание референсного генома...")
        urllib.request.urlretrieve("https://ctlab.itmo.ru/~aivanov/ct_algo_2026/ecoli_1k.fna", ref_file)
    if not os.path.exists(reads_file):
        print("Скачивание прочтений...")
        urllib.request.urlretrieve("https://ctlab.itmo.ru/~aivanov/ct_algo_2026/ecoli_reads.fastq", reads_file)
    os.makedirs(data_dir, exist_ok=True)

    k_values = [15, 21, 31]

    print("\n=== Обработка референсного генома ===")
    ref_seqs = read_fasta(ref_file)
    for k in k_values:
        print(f"Построение графа для k={k}...")
        graph = DeBruijnGraph(k)
        for seq in ref_seqs:
            graph.add_sequence(seq)
        graph.build_structure()
        graph.compress()
        graph.save_contigs(os.path.join(data_dir, f"ref_contigs_k{k}.fasta"))
        graph.save_gfa(os.path.join(data_dir, f"ref_graph_k{k}.gfa"))
        print(f"  Статистика: {graph.get_stats()}")

    print("\n=== Обработка прочтений ===")
    read_seqs = read_fastq(reads_file)
    for k in k_values:
        print(f"Построение графа для k={k}...")
        graph = DeBruijnGraph(k)
        for seq in read_seqs:
            graph.add_sequence(seq)

        graph.build_structure()
        print(f"  До очистки: {graph.get_stats()}")
        graph.compress()
        graph.save_contigs(os.path.join(data_dir, f"reads_uncleaned_contigs_k{k}.fasta"))
        graph.save_gfa(os.path.join(data_dir, f"reads_uncleaned_graph_k{k}.gfa"))

        print(f"  Очистка и повторное сжатие...")
        graph.clean_uncompressed(min_cov=1, max_tip_len_k1=5)
        print(f"  После очистки: {graph.get_stats()}")
        graph.save_contigs(os.path.join(data_dir, f"reads_cleaned_contigs_k{k}.fasta"))
        graph.save_gfa(os.path.join(data_dir, f"reads_cleaned_graph_k{k}.gfa"))

    print("\n\n=== Отрисовка графов ===")
    if not os.path.exists(bandage_exc):
        print("Отрисовать графы невозможно, так как передане некорректный путь к программе Bandage")
    else:
        for file in os.listdir(data_dir):
            if not file.endswith(".gfa"):
                continue
            print(f"Обработка файла {file:30}...", end="")
            subprocess.run([
                bandage_exc, "image",
                os.path.join(data_dir, file),
                os.path.join(data_dir, file.replace(".gfa", ".png"))
            ], capture_output=True, text=True)
            print(" готво!")


if __name__ == "__main__":
    main()