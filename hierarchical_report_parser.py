
import csv
import itertools

class ParsedRow(object):

    row = None

    def get_line_number(self):
        return self.row[22]

    def __str__(self):

        return "{0}(line={1})".format(self.__class__.__name__, self.get_line_number())

class ProteinRow(ParsedRow):


    peptide_rows = None
    accession = None
    validation = None
    confidence_of_specificity = None


    mw = None
    possible_coverage_percent = None
    coverage_percent = None
    spectrum_counting_nsaf = None
    confidence = None

    prot_id_lookup_map = {}


    def related_protein_seq_ids(self, is_legacy):

        res = []

        for p_accession in self.related_proteins:
            prot_id = self.prot_id_lookup_map.get(p_accession)

            if prot_id is not None:
                res.append(prot_id)
            elif not is_legacy:
                print("unknown prot: {0}".format(p_accession))

        return list(set(res))

    keep = False
    single_gene_origin_to_validate = False

    def __init__(self, row, peptide_rows):

        self.row = row
        self.peptide_rows = peptide_rows
        self.accession = row[1]
        self.validation = row[21]
        self.mw = float(row[3])
        self.possible_coverage_percent = float(row[4])
        self.coverage_percent = float(row[5])
        self.spectrum_counting_nsaf = float(row[6])
        self.confidence = float(row[20])

        self.related_proteins = list(map(lambda s: s.strip(), row[13].split(",")))

class PeptideRow(ParsedRow):

    psm_rows = None
    peptide_names = None
    keep = False

    def __init__(self, row, psm_rows):

        self.row = row
        self.psm_rows = psm_rows
        self.peptide_names = row[1].split(";")

"""

    parses :

     "Oxidation of M (6: Very Confident, 9: Very Confident) Oxidation of M (4: 100.0) Lysine 13C( 6) 15N(2) (13: Very Confident)"

    into :

     ['Oxidation of M (6: Very Confident, 9: Very Confident)',
      'Oxidation of M (4: 100.0)',
      'Lysine 13C( 6) 15N(2) (13: Very Confident)'
     ]


     "Phosphorylation of S (Not Scored), Phosphorylation of T (3: Doubtfull), Phosphorylation of Y (Not Scored)"

     into :

     ['Phosphorylation of T (3: Doubtfull']


"""
def split_mods(s):

    state = 0
    idx = 0
    p0 = 0
    tot = len(s)

    for c in s:

        if state == 0 and c == "(":
            if s[idx:idx+12] == "(Not Scored)":
                if idx+12 == tot:
                    # the last mod is Not Scored, get out...
                    break
                p0 = idx+ 13

        if state == 0 and c == ":":
            state = 1

        if state == 1 and c == ")":

            state = 0
            res = s[p0:idx]
            p0 = idx+1
            yield res.strip()

        idx += 1

class PSMRow(ParsedRow):

    protein_names = None
    validation = None
    variable_modifications = None
    modified_sequence = None
    sequence = None
    confidence  = None
    mz = None
    spectrum_scan_number = None
    spectrum_title = None
    retention_time = None

    single_gene_origin_to_validate = False

    confidence_of_specificity = None

    modifications = []

    confident_modifications = []

    keep = False

    def __init__(self, row):

        self.row = row
        self.protein_names = map(lambda n: n.strip(), row[1].split(","))
        self.sequence = row[2]
        self.modified_sequence = row[3]
        self.validation = row[21]
        self.variable_modifications = row[4]
        self.confidence  = float(row[20])
        self.mz = float(row[10])
        self.precursor_mz_error = float(row[15])

        self.retention_time = float(row[9])

        #if '=' in row[7]:
        #    self.spectrum_scan_number = int(row[7].split("=")[1])

        self.spectrum_file = row[6]
        self.spectrum_title = row[7]


        def parse_S_col(col_idx):

            col = row[col_idx]

            #print("col{0}: {1}".format(col_idx, col))

            for m0 in split_mods(col):
                if m0 == '':
                    continue

                mod_type_0, mods = m0.rsplit("(",1)
                mod_type = mod_type_0.strip(" ,")

                #if len(mods.split(",")) > 1:
                #    print("----> {0}".format(self.get_line_number()))

                for mod in mods.split(","):

                    if mod == "Not Scored":
                        continue

                    try:
                        pos, conf_score = mod.split(":")

                        pos = int(pos)

                        yield [col_idx, mod_type, pos, conf_score.strip()]
                    except ValueError as ve:
                        raise ve


        def merge_confidence_and_confidence_score_from_col_17_and_18():

            condifence_by_pos = dict([(m[2], m) for m in parse_S_col(17)])

            condifence_score_by_pos = [(m[2], m) for m in parse_S_col(18)]

            for pos, col_18 in condifence_score_by_pos:

                confidence = condifence_by_pos.get(pos)

                res = {
                    "modification_type": col_18[1],
                    "position": col_18[2],
                    "confidence": None if confidence is None else confidence[3],
                    "confidence_score": float(col_18[3])
                }

                yield res



        self.modifications = list(merge_confidence_and_confidence_score_from_col_17_and_18())




def __array_is_strict_prefix_of(a, b):

    if len(a) >= len(b) or len(a) == 0:
        return False

    return a == [
        a0 for a0, b0 in zip(a, b) if a0 == b0
    ]




def test_array_is_strict_prefix_of():

    def gen():
        yield [1,2], [1,2,1], True
        yield [1,2], [1,1,1], False
        yield [1,2], [2,2,1], False
        yield [1,2], [1,2], False
        yield [1,2], [2,2], False
        yield [2,2], [2,2,4], True
        yield [5], [5,2,4], True
        yield [2,3,5], [2], False
        yield [5], [5], False
        yield [], [], False
        yield [], [2], False
        yield [1], [], False
        yield [1,2,1], [1,1,1,1], False

    for a1, a2, expected in gen():
        print(a1, a2, expected)
        r = __array_is_strict_prefix_of(a1, a2)

        if r != expected:
            raise Exception("expected {0} for array_is_strict_prefix_of({1}, {2}), got {3}".format(expected, a1, a2, r))





def parse(complete_file_name, stop_at_line_number=None):

    """
    :param complete_file_name:
    :param stop_at_line_number: if None, parses the whole file, if not None, stops parsing line "stop_at_line_number"
    :return: an iterator of ProteinRow
    """

    c = 0

    def tsv_iterator():
        with open(complete_file_name) as file:
            tsvin = csv.reader(file, delimiter='\t')
            c = 0
            for row in tsvin:
                if stop_at_line_number is not None and c >= stop_at_line_number:
                    break

                c += 1


                if row[0] == '':
                    continue

                # replace string "path" with list (higher level) :
                row[0] = row[0].split(".")

                # add line number  at the end for debugging :
                row.append(c)

                yield row

    for _, tree_rows in itertools.groupby(tsv_iterator(), lambda row: row[0][0]):

        # root is always first
        root_row = next(tree_rows)

        # children follow
        children_rows = list(tree_rows)

        peptide_rows = []
        psm_rows = []

        for child_row in children_rows:
            if len(child_row[0]) == 2:
                peptide_rows.append(child_row)
            elif len(child_row[0]) == 3:
                psm_rows.append(child_row)


        yield ProteinRow(root_row, [
            PeptideRow(
                peptide_row,
                [ PSMRow(psm)
                  for psm in filter(lambda r: __array_is_strict_prefix_of(peptide_row[0], r[0]), psm_rows)
                ]
            )
            for peptide_row in peptide_rows
        ])


def test_parse2():

    c = 0

    for p in parse("/vagrant/tmp/ANTONV_J130605_005_sample_1_Default_Hierarchical_Report.tsv"):
        print("prot ", c)
        c += 1

def test_parse3(f):
    for p in parse(f):
        for pep in p.peptide_rows:
            __array_is_strict_prefix_of(p.row[0], pep.row[0])
            for psm in pep.psm_rows:
                __array_is_strict_prefix_of(pep.row[0], psm.row[0])


#test_parse4()
#test_array_is_strict_prefix_of()

#test_parse3("/vagrant/GM18486_MSB15024_01B_sample_1_Default_Hierarchical_Report.txt")


#print(list(split_mods("Oxidation of M (6: Very Confident, 9: Very Confident) Oxidation of M (4: 100.0) Lysine 13C( 6) 15N(2) (13: Very Confident)")))

#print(str(list(split_mods(
#    "Oxidation of M (6: Very Confident, 9: Very Confident) Oxidation of M (4: 100.0) Lysine 13C( 6) 15N(2) (13: Very Confident)"
#))))

#print(str(list(split_mods(
#    "Phosphorylation of S (Not Scored), Phosphorylation of T (3: Doubtfull), Phosphorylation of Y (Not Scored) Oxidation of M (6: Very Confident, 9: Very Confident) Oxidation of M (4: 100.0) Lysine 13C( 6) 15N(2) (13: Very Confident)"
#))))
