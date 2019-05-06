import mysql.connector
import pybedtools
import pysam

__author__ = 'Stefan Brenner'

class Assignment1:
    
    def __init__(self):
        ## Your gene of interest
        self.gene = "HLCS"
        self.gene_dictionary = self.download_gene_coordinates(file_name="UCSC_database.txt")
        self.bamfile = "chr21.bam"
        self.samfile = pysam.AlignmentFile(self.bamfile, "rb")
        self.reads = list(self.samfile.fetch("chr21", self.gene_dictionary["txStart"], self.gene_dictionary["txEnd"]))

    
    def download_gene_coordinates(self, file_name):
        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db="hg38")
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields) + ' WHERE refGene.name2="' + self.gene + '"'
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        gene_dictionary = {}
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
                gene_dictionary = {
                    "name2": row[0],
                    "name": row[1],
                    "chrom": row[2],
                    "txStart": row[3],
                    "txEnd": row[4],
                    "strand": row[5],
                    "exonCount": row[6],
                    "exonStarts": row[7],
                    "exonEnds": row[8]
                }
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        print("Done fetching data")
        return gene_dictionary
        
    def get_coordinates_of_gene(self):
        print("Coordinates_of_gene " + self.gene)
        print("     Start:\t", self.gene_dictionary["txStart"], "\n     End:\t", self.gene_dictionary["txEnd"],"\n")
        
    def get_gene_symbol(self):
        print("Gene Symbol: ",self.gene_dictionary["name"],"\n")
                        
    def get_sam_header(self):
        print("Samfile Header","\n", self.samfile.header['HD'])

    def get_properly_paired_reads_of_gene(self):
        counting_reads = 0
        for found_read in self.reads:
            if found_read.is_proper_pair:
                counting_reads += 1
        print("Properly paired reads: ",counting_reads,"\n")

    def get_gene_reads_with_indels(self):
        counting_reads = 0
        for found_read in self.reads:
            if not found_read.is_unmapped:
                cigar_search = found_read.cigar
                for (type, length) in cigar_search:
                    if (type == 1) or (type == 2):
                        counting_reads += 1
        print("Gene reads with indels:",counting_reads,"\n")
        
    def calculate_total_average_coverage(self):
        open_with_bedtools = pybedtools.BedTool(self.bamfile)
        genome_coverage = open_with_bedtools.genome_coverage(bg=True)

        total = 0
        average = 0

        for line in genome_coverage:
            number = float(line[3])
            average += number
            total += 1

        coverage = average / total
        print("Total average coverage:", coverage,"\n")
        
    def calculate_gene_average_coverage(self):
        open_with_bedtools = pybedtools.BedTool(self.bamfile)
        genome_coverage = open_with_bedtools.genome_coverage(bg=True)

        total = 0
        average = 0

        for line in genome_coverage:
            number = float(line[3])
            cbeg = int(line[1])

            if cbeg > self.gene_dictionary["txStart"]:
                if int(line[2]) <= self.gene_dictionary["txEnd"]:
                    average += number
                    total += 1

        coverage = average / total
        print("Total gene average coverage:", coverage,"\n")
        
    def get_number_mapped_reads(self):
        counting_reads = 0
        for found_read in self.reads:
            if not found_read.is_unmapped:
                counting_reads += 1
        print("Mapped reads:",counting_reads,"\n")

    def get_region_of_gene(self):
        print("Region of gene: \nChromosome ", self.gene_dictionary["chrom"],"\n")

    def get_number_of_exons(self):
        print("Number of exons: ", self.gene_dictionary["exonCount"],"\n")
    
    def print_summary(self):
        self.get_coordinates_of_gene()
        self.get_gene_symbol()
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_total_average_coverage()
        self.calculate_gene_average_coverage()
        self.get_number_mapped_reads()
        self.get_region_of_gene()
        self.get_number_of_exons()

    
def main():
    print("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()
    print("Done with assignment 1")

        
if __name__ == '__main__':
    main()
    
    
