#!/usr/bin/env python

"""
Script to convert the FH VCF to a CSV file and filter at a given gnomad pop frequency.

Usage - 

python fh_reporter.py --vcf $INPUT_VCF \
--sample_id $SAMPLE_ID \
--worklist_id $WORKLIST_ID \
--output $OUTPUT_LOCATION \
--gnomad_max_af $MAX_POP_FREQ

Requirements:

Input Files:

VCF - GATK multisample VCF normalised and with split multiallelics and annotated with VEP.

Software:

- python3
- pandas
- pyvariantfilter

TODO - Add quality filter for initial import!

"""


import argparse
import pandas as pd
from pyvariantfilter.family import Family
from pyvariantfilter.family_member import FamilyMember
from pyvariantfilter.variant_set import VariantSet

parser = argparse.ArgumentParser(description='Report FH variants in the custom ROIs')
parser.add_argument('--vcf', type=str, nargs=1, required=True,
				help='The input location of the multisample GATK VCF.')
parser.add_argument('--sample_id', type=str, nargs=1,
				help='The sample_id to report.')
parser.add_argument('--worklist_id', type=str, nargs=1,
				help='The worklist_id to include in the report.')
parser.add_argument('--output', type=str, nargs=1,
				help='The location and name of the output file.')
parser.add_argument('--gnomad_max_af', type=float, nargs=1,
				help='The maximium population frequency in gnomad for filtering on.')
args = parser.parse_args()

# Utility Functions

def passes_initial_filter(variant, proband_id):
	"""
	Filter variants from the VCF.
	
	We import if the variant passes quality filtering.
  
	"""
	if variant.has_alt(proband_id) and variant.passes_gt_filter(proband_id, min_gq=20) and variant.passes_filter():
		
		return True
		
	return False



def passes_final_filter(variant, proband_id):
	"""
	Filter variants from the VCF.
	
	We import if the variant passes quality filtering and is below 1% in gnomad exomes and gnomad genomes.
	
	"""

	# If the proband has the variant and we pass the genotype and variant level filters
	if variant.has_alt(proband_id) and variant.passes_gt_filter(proband_id, min_gq=20) and variant.passes_filter():
		
		
		# The filter_on_numerical_transcript_annotation_lte() function allows us to filter on numerical values 
		# we can set different cutoffs for different variant types. For example ad_het is variants in which the 
		# proband is heterozygous on an autosome. In this case we get two boolean values describing whether the 
		# variant is below x% in the gnomad genomes and gnomad exomes datasets.
		freq_filterg = variant.filter_on_numerical_transcript_annotation_lte(annotation_key='gnomADg_AF_POPMAX',
																						  ad_het=gnomad_max_af,
																						  ad_hom_alt=gnomad_max_af,
																						  x_male =gnomad_max_af,
																						  x_female_het=gnomad_max_af,
																						  x_female_hom=gnomad_max_af,
																						  compound_het=gnomad_max_af,
																						  y=gnomad_max_af,
																						  mt=gnomad_max_af,
																						  )
		freq_filtere = variant.filter_on_numerical_transcript_annotation_lte(annotation_key='gnomADe_AF_POPMAX',
																						  ad_het=gnomad_max_af,
																						  ad_hom_alt=gnomad_max_af,
																						  x_male =gnomad_max_af,
																						  x_female_het=gnomad_max_af,
																						  x_female_hom=gnomad_max_af,
																						  compound_het=gnomad_max_af,
																						  y=gnomad_max_af,
																						  mt=gnomad_max_af,
																						  )  
		
		
		if freq_filtere and freq_filterg:

			return True
		
	return False



# Create python objects for turning to dataframe.

sample_id = args.sample_id[0]
worklist_id = args.worklist_id[0]
vcf = args.vcf[0]
output = args.output[0]
gnomad_max_af = args.gnomad_max_af[0]

# Create a Family object - assume male as we don't have any targets on sex chromosomes so sex doesn't matter.
proband = FamilyMember(sample_id, 'FAM001', 1, True)
my_family = Family('FAM001')
my_family.add_family_member(proband)
my_family.set_proband(proband.get_id())

# Create a new VariantSet object
my_variant_set = VariantSet()

# Associate the my_family object with my_variant_set
my_variant_set.add_family(my_family)

# Read in variants from vcf
my_variant_set.read_variants_from_platypus_vcf(vcf,
 									proband_variants_only=True,
									filter_func=passes_initial_filter,
									args=(proband.get_id(),))

# Convert to dataframe
df_unfiltered = my_variant_set.to_df()

columns_to_keep = [
   '#SampleId',
   'WorklistId',
   'Variant',
   'Genotype',
   'WorstConsquence',
   'Feature',
   'ExistingVariation',
   'GnomadePopMax',
   'GnomadgPopMax',
   'Symbol',
   'Consequence',
   'HGVSp',
   'HGVSc',
   'Clinvar',
   'AutoTranscriptPick',
   'SIFT',
   'PolyPhen',
   'Exon',
   'Intron'

]

# Write emptry dataframe to csv if empty
if df_unfiltered.shape[1] == 0:

	df_unfiltered = pd.DataFrame(columns = columns_to_keep)

	df_unfiltered.to_csv(output +  '_unfiltered.csv', index=False, sep='\t')

else:

	# Rename columns and write to csv
	df_unfiltered['#SampleId'] = sample_id
	df_unfiltered['WorklistId'] = worklist_id
	df_unfiltered['Variant'] = df_unfiltered['variant_id']
	df_unfiltered['Genotype'] = df_unfiltered[f'{sample_id}_GT']
	df_unfiltered['WorstConsquence'] = df_unfiltered['worst_consequence']
	df_unfiltered['Feature'] = df_unfiltered['csq_Feature']
	df_unfiltered['ExistingVariation'] = df_unfiltered['csq_Existing_variation']
	df_unfiltered['GnomadePopMax'] = df_unfiltered['csq_gnomADe_AF_POPMAX']
	df_unfiltered['GnomadgPopMax'] = df_unfiltered['csq_gnomADg_AF_POPMAX']
	df_unfiltered['Symbol'] = df_unfiltered['csq_SYMBOL']
	df_unfiltered['Consequence'] = df_unfiltered['csq_Consequence']
	df_unfiltered['HGVSp'] = df_unfiltered['csq_HGVSp']
	df_unfiltered['HGVSc'] = df_unfiltered['csq_HGVSc']
	df_unfiltered['Clinvar'] = df_unfiltered['csq_CLIN_SIG']
	df_unfiltered['AutoTranscriptPick'] = df_unfiltered['csq_PICK']
	df_unfiltered['SIFT'] = df_unfiltered['csq_SIFT']
	df_unfiltered['PolyPhen'] = df_unfiltered['csq_PolyPhen']
	df_unfiltered['Exon'] = df_unfiltered['csq_EXON']
	df_unfiltered['Intron'] = df_unfiltered['csq_INTRON']

	df_unfiltered[columns_to_keep].to_csv(output +  '_unfiltered.csv', index=False, sep='\t')


# Now filter and write to csv
my_variant_set.filter_variants(passes_final_filter, args=(sample_id,))

df_filtered = my_variant_set.to_df()

# Write emptry dataframe to csv if empty
if df_filtered.shape[1] == 0:

	df_filtered = pd.DataFrame(columns = columns_to_keep)

	df_filtered.to_csv(output +  '_filtered.csv', index=False, sep='\t')

else:

	# Rename columns and write to csv
	df_filtered['#SampleId'] = sample_id
	df_filtered['WorklistId'] = worklist_id
	df_filtered['Variant'] = df_filtered['variant_id']
	df_filtered['Genotype'] = df_filtered[f'{sample_id}_GT']
	df_filtered['WorstConsquence'] = df_filtered['worst_consequence']
	df_filtered['Feature'] = df_filtered['csq_Feature']
	df_filtered['ExistingVariation'] = df_filtered['csq_Existing_variation']
	df_filtered['GnomadePopMax'] = df_filtered['csq_gnomADe_AF_POPMAX']
	df_filtered['GnomadgPopMax'] = df_filtered['csq_gnomADg_AF_POPMAX']
	df_filtered['Symbol'] = df_filtered['csq_SYMBOL']
	df_filtered['Consequence'] = df_filtered['csq_Consequence']
	df_filtered['HGVSp'] = df_filtered['csq_HGVSp']
	df_filtered['HGVSc'] = df_filtered['csq_HGVSc']
	df_filtered['Clinvar'] = df_filtered['csq_CLIN_SIG']
	df_filtered['AutoTranscriptPick'] = df_filtered['csq_PICK']
	df_filtered['SIFT'] = df_filtered['csq_SIFT']
	df_filtered['PolyPhen'] = df_filtered['csq_PolyPhen']
	df_filtered['Exon'] = df_filtered['csq_EXON']
	df_filtered['Intron'] = df_filtered['csq_INTRON']

	df_filtered[columns_to_keep].to_csv(output +  '_filtered.csv', index=False, sep='\t')