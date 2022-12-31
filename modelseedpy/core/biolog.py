# -*- coding: utf-8 -*-
pm1_str = """A01    Negative Control
A02    L-Arabinose    cpd00224
A03    N-Acetyl-DGlucosamine    cpd00122
A04    D-Saccharic Acid    cpd00571
A05    Succinic Acid    cpd00036
A06    D-Galactose    cpd00108
A07    L-Aspartic Acid    cpd00041
A08    L-Proline    cpd00129
A09    D-Alanine    cpd00117
A10    D-Trehalose    cpd00794
A11    D-Mannose    cpd00138
A12    Dulcitol     cpd01171
B01    D-Serine    cpd00550
B02    D-Sorbitol    cpd00588
B03    Glycerol    cpd00100
B04    L-Fucose    cpd00751
B05    D-Glucuronic Acid
B06    D-Gluconic Acid
B07    D,L-α-GlycerolPhosphate
B08    D-Xylose    cpd00154
B09    L-Lactic Acid    cpd00159
B10    Formic Acid    cpd00047
B11    D-Mannitol    cpd00314
B12    L-Glutamic Acid    cpd00023
C01    Glucose-6-Phosphate    cpd00079
C02    D-Galactonic Acid-γ-Lactone    cpd02143
C03    D,L-Malic Acid    cpd00130;cpd00386
C04    D-Ribose    cpd00105
C05    Tween 20    cpd13391;cpd01741
C06    L-Rhamnose    cpd00396
C07    D-Fructose    cpd00082
C08    Acetic Acid    cpd00029
C09    α-D-Glucose    cpd00027
C10    Maltose    cpd00179
C11    D-Melibiose    cpd03198
C12    Thymidine     cpd00184
D01    L-Asparagine    cpd00132
D02    D-Aspartic Acid    cpd00320
D03    D-Glucosaminic Acid    cpd02351
D04    1,2-Propanediol    cpd01861
D05    Tween 40    cpd24417;cpd00214
D06    α-Keto-Glutaric Acid    cpd00024
D07    α-Keto-Butyric Acid    cpd00094
D08    α-Methyl-DGalactoside
D09    α-D-Lactose    cpd00208
D10    Lactulose    cpd04349
D11    Sucrose    cpd00076
D12    Uridine     cpd00249
E01    L-Glutamine    cpd00053
E02    M-Tartaric Acid    cpd00432
E03    Glucose-1-Phosphate    cpd00089
E04    Fructose-6-Phosphate    cpd00072
E05    Tween 80    cpd13392
E06    α-Hydroxy Glutaric Acid-γLactone    cpd23862
E07    α-Hydroxy Butyric Acid
E08    β-Methyl-DGlucoside
E09    Adonitol    cpd00366
E10    Maltotriose    cpd01262
E11    2-Deoxy Adenosine    cpd00115
E12    Adenosine     cpd00182
F01    Glycyl-L-Aspartic Acid    cpd11589
F02    Citric Acid    cpd00137
F03    M-Inositol    cpd00121
F04    D-Threonine    cpd00611
F05    Fumaric Acid    cpd00106
F06    Bromo Succinic Acid    cpd23859
F07    Propionic Acid
F08    Mucic Acid
F09    Glycolic Acid    cpd00139
F10    Glyoxylic Acid    cpd00040
F11    D-Cellobiose
F12    Inosine     cpd00246
G01    Glycyl-L-Glutamic Acid
G02    Tricarballylic Acid    cpd16654
G03    L-Serine    cpd00054
G04    L-Threonine    cpd00161
G05    L-Alanine    cpd00035
G06    L-Alanyl-Glycine
G07    Acetoacetic Acid    cpd00142
G08    N-Acetyl-β-D-Mannosamine
G09    Mono Methyl Succinate
G10    Methyl Pyruvate    cpd24420
G11    D-Malic Acid    cpd00386
G12    L-Malic Acid    cpd00130
H01    Glycyl-L-Proline
H02    P-Hydroxy Phenyl Acetic Acid
H03    M-Hydroxy Phenyl Acetic Acid
H04    Tyramine    cpd00374
H05    D-Psicose    cpd03884
H06    L-Lyxose    cpd01067
H07    Glucuronamide    cpd26269
H08    Pyruvic Acid    cpd00020
H09    L-Galactonic Acid-γ-Lactone
H10    D-Galacturonic Acid
H11    Phenylethylamine    cpd03161
H12    2-Aminoethano    cpd00162"""


class BiologPlate:
    def __init__(self, plate_id, rows, cols):
        self.id = plate_id
        self.rows = rows
        self.cols = cols
        self.wells = {}
        self.base = {}

    def add_base(self, compounds, value):
        for c in compounds:
            self.base[c] = value

    def get_media(self, well_id):
        media = {}
        if well_id in self.wells:
            well = self.wells[well_id]
            for compound in self.base:
                media[compound] = self.base[compound]
            for compound in well["compounds"]:
                media[compound] = well["value"]
            return media
        return None

    def _repr_html_(self):
        html = "<h3>" + self.id + "</h3>"
        t = '<table style="border-width: 1px; box-sizing: border-box; border-color: grey"><tbody>'
        for row in self.rows:
            t += "<tr>"
            for col in self.cols:
                t += '<td style="text-align:left;vertical-align:top;">'
                for o in self.wells:
                    if o == row + col:
                        t += self.wells[o]["desc"]
                        # print(self.wells[o])
                t += "</td>"
            t += "</tr>"
        t += "</tbody></table>"
        html += t
        return html


class Biolog:
    def __init__(self):
        self.plates = {}

    def add_plate(self, plate):
        if plate.id in self.plates:
            print("replace existing plate")
        self.plates[plate.id] = plate

    def run_plates(self, model, biomass=None, cmp="e"):  # !!! biomass is never used
        prev_medium = model.medium
        compound_exchange = {}
        for ex_rxn in model.exchanges:
            ex_met = list(ex_rxn.metabolites)[0]
            if ex_met.id not in compound_exchange:
                compound_exchange[ex_met.id] = ex_rxn
            else:
                print(ex_rxn, list(ex_rxn.metabolites)[0])
        for plate_id, plate in self.plates.items():
            for well_id in plate.wells:
                media = plate.get_media(well_id)
                model_medium = {}
                for o in media:
                    match = "{}_{}".format(o, cmp)
                    if match in compound_exchange:
                        model_medium[compound_exchange[match].id] = media[o]
                    else:
                        print("skip", o)
                model.medium = model_medium
                ovalue = model.slim_optimize()
                plate.wells[well_id]["growth"] = ovalue
                # print(well_id, solution)
                # print(well_id, media)
        return prev_medium
