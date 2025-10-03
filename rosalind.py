from collections import OrderedDict
from enum import Enum
import itertools

def read_rosalind(name: str) -> list[str]:
    with open(f"rosalind_{name}.txt") as f:
        lines = f.readlines()
        return [line.strip() for line in lines]

COMPLEMENTS = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
}

RNA_CODON = {
    "UUU":  "F",
    "CUU":  "L",
    "AUU":  "I",
    "GUU":  "V",
    "UUC":  "F",
    "CUC":  "L",
    "AUC":  "I",
    "GUC":  "V",
    "UUA":  "L",
    "CUA":  "L",
    "AUA":  "I",
    "GUA":  "V",
    "UUG":  "L",
    "CUG":  "L",
    "AUG":  "M",
    "GUG":  "V",
    "UCU":  "S",
    "CCU":  "P",
    "ACU":  "T",
    "GCU":  "A",
    "UCC":  "S",
    "CCC":  "P",
    "ACC":  "T",
    "GCC":  "A",
    "UCA":  "S",
    "CCA":  "P",
    "ACA":  "T",
    "GCA":  "A",
    "UCG":  "S",
    "CCG":  "P",
    "ACG":  "T",
    "GCG":  "A",
    "UAU":  "Y",
    "CAU":  "H",
    "AAU":  "N",
    "GAU":  "D",
    "UAC":  "Y",
    "CAC":  "H",
    "AAC":  "N",
    "GAC":  "D",
    "UAA":  "Stop",
    "CAA":  "Q",
    "AAA":  "K",
    "GAA":  "E",
    "UAG":  "Stop",
    "CAG":  "Q",
    "AAG":  "K",
    "GAG":  "E",
    "UGU":  "C",
    "CGU":  "R",
    "AGU":  "S",
    "GGU":  "G",
    "UGC":  "C",
    "CGC":  "R",
    "AGC":  "S",
    "GGC":  "G",
    "UGA":  "Stop",
    "CGA":  "R",
    "AGA":  "R",
    "GGA":  "G",
    "UGG":  "W",
    "CGG":  "R",
    "AGG":  "R",
    "GGG":  "G",
}

class SequenceType(Enum):
    UNKNOWN = {}
    DNA = {'A', 'C', 'G', 'T'}
    RNA = {'A', 'C', 'G', 'U'}
    PROTEIN = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}

class Sequence(str):
    __sequence_type: SequenceType = SequenceType.UNKNOWN

    def __new__(cls, value: str | None, sequence_type: SequenceType = SequenceType.UNKNOWN):
        obj = super().__new__(cls, value)
        obj.__sequence_type = sequence_type if sequence_type != SequenceType.UNKNOWN else obj.__infer_sequence_type()

        return obj

    @property
    def sequence_type(self) -> SequenceType:
        return self.__sequence_type

    @property
    def symbols_count(self) -> dict[str, int]:
        counts = {}

        for char in self:
            counts[char] = counts.get(char, 0) + 1

        return counts
    
    @property
    def gc_ratio(self) -> float:
        counts = self.symbols_count
        percentage = (counts['C'] + counts['G']) / len(self) * 100
        return percentage
    
    @property
    def complement(self) -> 'Sequence':
        result = ""
        for char in self:
            result = COMPLEMENTS[char] + result
        
        return Sequence(result)
    
    @property
    def protein_string(self) -> 'Sequence':
        if self.sequence_type == SequenceType.RNA:
            codon_table = RNA_CODON
        elif self.sequence_type == SequenceType.DNA:
            raise NotImplementedError("protein_string not implemented for DNA sequence")
        else:
            raise ValueError(f"Unable to convert sequence of type {self.__sequence_type.name} to protein")
        
        result = Sequence("", sequence_type=SequenceType.PROTEIN)
        for i in range(len(self) // 3):
            group = self[i*3:(i+1)*3]
            symbol = codon_table[group]

            if symbol == 'Stop':
                break

            result += symbol

        return Sequence(result)
    
    def sub_locations(self, other: 'Sequence') -> list[int]:
        result = []
        for i in range(0, len(self) - len(other)):
            if self[i:i+len(other)] == other:
                result.append(i+1)

        return result

    def hamming_distance(self, other: 'Sequence') -> int:
        if len(self) != len(other):
            raise ValueError(f"Length of sequences differ. Self: {len(self)}, other: {len(other)}")

        result = 0

        for i in range(len(self)):
            result += 1 if self[i] != other[i] else 0

        return result

    def overlaps(self, other: 'Sequence', length: int) -> bool:
        return self[-length:] == other[:length]
    
    def __infer_sequence_type(self) -> SequenceType:
        symbols = self.symbols_count.keys()

        for sequence_type in SequenceType:
            if symbols - sequence_type.value == set():
                return sequence_type
            
        return SequenceType.UNKNOWN

class FASTA(OrderedDict[str | None, Sequence]):
    @property
    def sequence(self) -> Sequence:
        return self[None]
    
    def overlap_graph(self, length: int) -> list[tuple[str, str]]:
        result = []

        for (id_a, sequence_a), (id_b, sequence_b) in itertools.product(self.items(), self.items()):
            if sequence_a == sequence_b:
                continue

            if sequence_a.overlaps(sequence_b, length):
                result.append((id_a, id_b))

        return result


    @staticmethod
    def parse(text) -> 'FASTA':
        result = OrderedDict()

        current_id = None
        current_sequence = ""
        reading_id = False

        for char in text:
            if char == ">":
                if current_id:
                    result[current_id] = Sequence(current_sequence)

                current_id = ""
                current_sequence = ""
                reading_id = True

                continue

            if char == "\n":
                reading_id = False
                continue

            if reading_id:
                current_id += char
                continue

            current_sequence += char

        result[current_id] = Sequence(current_sequence)

        return FASTA(result)

    @staticmethod
    def from_file(file_path: str) -> 'FASTA':
        with open(file_path) as f:
            return FASTA.parse(f.read())
        
    @staticmethod
    def from_rosalind(name: str) -> 'FASTA':
        return FASTA.from_file(f"rosalind_{name}.txt")
