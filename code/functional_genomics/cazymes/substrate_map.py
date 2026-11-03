"""Curated primary-substrate assignment for CAZy families (degradative classes).
Assignments follow the predominant characterised activity of each family; many families
are polyspecific, and families without a dominant substrate are 'Mixed / other'.
GTs are biosynthetic and kept as a separate class, not a carbon-utilisation substrate."""
import re

SUBSTRATE = {
 "Cellulose": "GH6 GH7 GH9 GH12 GH44 GH45 GH48 GH131 AA9 AA16 CBM1 CBM2 CBM3 CBM10 CBM63 CBM64",
 "Hemicellulose (xylan, xyloglucan, mannan, arabinan)":
   "GH10 GH11 GH26 GH39 GH43 GH51 GH54 GH62 GH67 GH74 GH113 GH115 GH120 GH127 GH134 GH142 GH146 "
   "CE1 CE2 CE3 CE5 CE6 CE7 CE15 CE16 AA14 CBM6 CBM22 CBM35 CBM42 CBM60 CBM29 CBM9 CBM4 CBM15",
 "Pectin": "GH28 GH53 GH78 GH88 GH93 GH105 GH106 GH138 GH140 PL1 PL3 PL4 PL9 PL10 PL11 PL22 PL26 CE8 CE12 CE13 CBM77",
 "Starch / glycogen": "GH13 GH14 GH15 GH57 GH77 GH119 GH126 GH133 AA13 CBM20 CBM21 CBM25 CBM26 CBM41 CBM48 CBM53 CBM69",
 "β-glucan / fungal & algal cell wall": "GH16 GH17 GH55 GH64 GH71 GH72 GH81 GH128 GH132 GH152 GH158 CBM43 CBM52 CBM39",
 "Chitin / chitosan": "GH18 GH19 GH46 GH75 CE4 CE9 AA10 AA11 AA15 CBM5 CBM12 CBM14 CBM18 CBM19 CBM54 CBM73",
 "Peptidoglycan (bacterial)": "GH23 GH24 GH25 GH73 GH103 GH104 GH108 CBM50",
 "Host / animal glycans (mucin, GAG, N-glycan)":
   "GH20 GH29 GH33 GH56 GH79 GH84 GH85 GH89 GH95 GH98 GH101 GH109 GH110 GH112 GH123 GH129 GH136 GH141 GH151 GH156 "
   "GH38 GH47 GH63 GH92 GH99 GH125 PL8 PL12 PL13 PL15 PL21 PL29 PL30 PL33 CE20 CBM32 CBM40 CBM47 CBM51 CBM67",
 "Sucrose / fructan / trehalose": "GH32 GH68 GH91 GH37 GH65 GH100",
 "α-galactosides / glycolipids": "GH27 GH36 GH4 GH59 GH97",
 "β-glucosides / β-galactosides (mixed)": "GH1 GH2 GH3 GH5 GH35 GH42 GH8 GH116 GH76 GH114 GH135 GH31 GH130 GH94",
 "Algal polysaccharides (alginate, carrageenan, agar, ulvan)": "PL6 PL7 PL14 PL17 PL18 PL24 PL25 PL28 PL36 PL38 PL40 GH82 GH86 GH117 GH118 GH150 GH167",
 "Lignin / oxidative (AA1–AA8, AA12)": "AA1 AA2 AA3 AA4 AA5 AA6 AA7 AA8 AA12",
 "Cellulosome (dockerin/cohesin)": "fungi_doc fungi_coh dockerin cohesin",
}
FAM2SUB = {f: s for s, fl in SUBSTRATE.items() for f in fl.split()}

def family(model: str) -> str:
    if model.startswith("GT2_"):
        return "GT2"
    m = re.match(r"^((?:GH|GT|PL|CE|AA|CBM)\d+)(?:_\d+)?$", model)
    return m.group(1) if m else model

def cazy_class(fam: str) -> str:
    m = re.match(r"^(GH|GT|PL|CE|AA|CBM)", fam)
    return m.group(1) if m else ("Transporter" if fam == "Sugar_tr" else "Other module")

def substrate(fam: str) -> str:
    if fam == "Sugar_tr":
        return "Sugar transporters (MFS Sugar_tr)"
    if fam.startswith("GT"):
        return "Glycosyltransferases (biosynthesis)"
    return FAM2SUB.get(fam, "Mixed / other")
