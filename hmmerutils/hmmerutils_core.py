import re
import os
import sys
import argparse
from multiprocessing import Pool
from hmmerutils.utils import fas_to_dic

# ref_genome.keys()
class Hmmer:

    def __init__(self, 
                 references = None,
                 tables = None, 
                 min_length = 100,
                 species = None,
                 threads = 1):
        
        self.references = references
        self.tables = tables
        self.min_length = min_length
        self.species = species
        self.threads = threads

        # placeholder (sequence in dict)
        self.scaffold = {}

    def revcom(self, string: str = None) -> str:
        """
        Reverse complement
        """
        libco = {"A": "T", "G": "C", "C": "G",
                "T": "A", "R": "Y", "Y": "R",
                "S": "S", "W": "W", "K": "M",
                "M": "K", "B": "V", "V": "B",
                "D": "H", "H": "D", "N": "N",
                "-": "-"}

        complement = []
        for pos in range(0, len(string)):
            complement.append(libco[string[pos]])

        return "".join(complement)[::-1]

    def select_slice(self, scaffold:str=None,positions:tuple=None,isrevcom:bool=False) -> str:

        pos1, pos2 = positions

        if isrevcom:
            return self.revcom( scaffold[pos2 - 1: pos1] )

        else:
            return scaffold[pos1 - 1: pos2]

    def gethits(self, table: str = None) -> dict:
        # table = tables[0]
        tarfilpatt = "# Target file:[ ]+?(.+)\n$"
        ref_file   = ''
        matches    = {}

        with open(table , 'r') as f:

            for line in f.readlines():

                if not line.startswith("#") and not matches:

                    colvals    = re.sub("[ ]+", " ", line).strip().split(' ')
                    pos1, pos2 = int(colvals[6]), int(colvals[7]) 

                    if abs(pos1 - pos2) + 1 < self.min_length:
                        continue

                    matches['locus']     = colvals[2]
                    matches['id']        = colvals[0]
                    matches['positions'] = (pos1, pos2)
                    matches['isrevcom']  = pos2 < pos1
                    matches['eval']      = float(colvals[12])
                    matches['score']     = float(colvals[13])

                else:

                    if line.startswith("# Target file:"):
                        ref_file += os.path.basename( re.sub(tarfilpatt, '\\1', line) )

        if matches:
            matches['ref_file'] = ref_file
            return matches

        else:
            return None

    def unfair_swap(self, reduced):

        preout = {}

        for k,v in reduced.items():
            ref_file = v['ref_file']

            if not preout.__contains__(ref_file):
                preout[ref_file]  = [ {k:v} ]

            else:
                preout[ref_file] += [ {k:v} ]

        return preout

    def reduce_hits(self, jsons):
        out = {}

        for json in jsons:
            locus = json['locus']
            load  = {k:v for k,v in json.items() if k != 'locus'}

            if not out.__contains__(locus):
                out[locus] = load

            else:
                if json['eval'] < out[locus]['eval']:
                    out[locus] = load

                elif json['eval'] == out[locus]['eval']:

                    if json['score'] > out[locus]['score']:
                        out[locus] = load

        return self.unfair_swap(out)

    def read_ref(self, reference):
        # reference = 'mock3.fna'
        complete_file = ''

        for i in self.references:
            basename = os.path.basename(i)
            if basename == reference:
                complete_file += i
                break
        
        if not complete_file:
            return None
        
        return fas_to_dic(complete_file)

    def export_exons(self, exon):

        exon_key  = list(exon.keys())[0]
        positions = exon[exon_key]['positions']
        isrevcom  = exon[exon_key]['isrevcom']

        idpatt = re.compile("^>%s " % exon[exon_key]['id'])
        #    output = exon_key
        for ref_id,ref_seq in self.scaffold.items():

            if idpatt.match(ref_id):
                tmp_seq = self.select_slice(ref_seq, positions, isrevcom)

                with open(exon_key, 'w') as f:
                    f.write( ">%s\n%s\n" % (self.species, tmp_seq) )

    def run(self):

        with Pool(processes = self.threads) as p:
        
            hitsmeta = [*p.map(self.gethits, self.tables)]
            jsons    = list(filter(None, hitsmeta))

            if not jsons:
                sys.stderr.write('No hit found\n')
                sys.stderr.flush()
                exit()

            reduced  = self.reduce_hits(jsons)

            for ref, exons in reduced.items():
                # ref = 'mock3.fna'
                sys.stderr.write('Slicing %s\n' % ref)
                sys.stderr.flush()

                self.scaffold = self.read_ref(ref)

                if not self.scaffold:
                    sys.stderr.write('"%s" doesn\'t match with any scaffold name\n' % ref)
                    sys.stderr.flush()
                    continue

                [*p.map(self.export_exons, exons)]

def getOpts():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
                           Parse HMMER output
    Example:

        $ hmmerparse -t [hit tables] -f [scaffolds]  -s [species mame]
""")

    parser.add_argument('-t','--tables',
                        required=True,
                        nargs='+',
                        metavar="",
                        action='store',
                        default=None,
                        help='Hit tables from HMMER')
    parser.add_argument('-f','--scaffolds',
                        required=True,
                        nargs='+',
                        metavar="",
                        action='store',
                        default=None,
                        help='Scaffolds or genomes sequences previously used with HMMER')
    parser.add_argument('-s','--species',
                        required=True,
                        metavar="",
                        type= str,
                        action='store',
                        default=None,
                        help='Species name or label to name sequences')
    parser.add_argument('-m','--min', 
                        metavar="",
                        type = int,
                        default = 100,
                        help='''[Optional] Minimum length for each hit [Default: %s] ''' % 100)
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default: 1]')
    args = parser.parse_args()
    return args

def main():
    args = getOpts()
    Hmmer(
        references = args.scaffolds,
        tables     = args.tables,
        min_length = args.min, 
        species    = args.species,
        threads    = args.threads
    ).run()

if __name__ == "__main__":
    main()

