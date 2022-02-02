import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("[OD Sample File]", help='tab delimited file with sampleID, dilution factor, and OD info')
args = parser.parse_args()

SampleFile = pd.read_csv(sys.argv[1], sep='\t')
print(type(SampleFile))
Double = []
for index, row in SampleFile.iterrows():
    DilutionFactor = row['DilutionFactor']
    TransferCultureOD = row['TransferCultureOD']
    InoculationOD: float = TransferCultureOD / DilutionFactor
    FinalOD = row['FinalOD']
    x = InoculationOD
    GenerationCounter = 0
    while x < FinalOD:
        GenerationCounter += 1
        x = x * 2
    LowGenerationOD = x / 2
    Double.append("%.3f" % float(((x - FinalOD) / LowGenerationOD) + (GenerationCounter - 1)))
SampleFile['doublings'] = Double
SampleFile.to_excel('AmpliconGenerations.xlsx')

print(SampleFile)


