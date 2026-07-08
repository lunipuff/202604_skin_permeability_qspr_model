############################################################
# calculate_rdkit_descriptors.py
# Calculate RDKit descriptors from curated SMILES table
############################################################

import os
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors


############################################################
# File paths
############################################################

input_path = "data/interim/compound_structure_lookup.csv"

output_as_reported_path = "data/processed/rdkit_descriptors_as_reported.csv"
output_parent_path = "data/processed/rdkit_descriptors_parent.csv"

log_path = "data/interim/rdkit_descriptor_calculation_log.csv"


############################################################
# Helper functions
############################################################

def safe_float(value):
	if value is None:
		return None

	try:
		return float(value)
	except Exception:
		return None


def calculate_descriptors(smiles_string):
	mol = Chem.MolFromSmiles(smiles_string)

	if mol is None:
		return None

	canonical_smiles = Chem.MolToSmiles(
		mol,
		isomericSmiles=True
	)

	descriptor_values = {
		"canonical_smiles": canonical_smiles,
		"rdkit_MolWt": safe_float(Descriptors.MolWt(mol)),
		"rdkit_ExactMolWt": safe_float(Descriptors.ExactMolWt(mol)),
		"rdkit_MolLogP": safe_float(Crippen.MolLogP(mol)),
		"rdkit_MolMR": safe_float(Crippen.MolMR(mol)),
		"rdkit_TPSA": safe_float(rdMolDescriptors.CalcTPSA(mol)),
		"rdkit_NumHDonors": safe_float(Lipinski.NumHDonors(mol)),
		"rdkit_NumHAcceptors": safe_float(Lipinski.NumHAcceptors(mol)),
		"rdkit_NumRotatableBonds": safe_float(Lipinski.NumRotatableBonds(mol)),
		"rdkit_RingCount": safe_float(Lipinski.RingCount(mol)),
		"rdkit_HeavyAtomCount": safe_float(Lipinski.HeavyAtomCount(mol)),
		"rdkit_FractionCSP3": safe_float(rdMolDescriptors.CalcFractionCSP3(mol)),
		"rdkit_NumAromaticRings": safe_float(rdMolDescriptors.CalcNumAromaticRings(mol)),
		"rdkit_NumAliphaticRings": safe_float(rdMolDescriptors.CalcNumAliphaticRings(mol)),
		"rdkit_NumSaturatedRings": safe_float(rdMolDescriptors.CalcNumSaturatedRings(mol)),
		"rdkit_FormalCharge": safe_float(Chem.GetFormalCharge(mol))
	}

	return descriptor_values


def calculate_descriptor_table(df, smiles_col, representation_name):
	results = []
	log_rows = []

	for _, row in df.iterrows():
		compound_id = row["compound_id"]
		cas_no = row["CAS.No"]
		compound = row["Compound"]
		smiles_string = row[smiles_col]

		base_log = {
			"compound_id": compound_id,
			"CAS.No": cas_no,
			"Compound": compound,
			"representation": representation_name,
			"smiles_column": smiles_col,
			"SMILES": smiles_string
		}

		if pd.isna(smiles_string) or str(smiles_string).strip() == "":
			log_row = dict(base_log)
			log_row["rdkit_status"] = "failed"
			log_row["rdkit_note"] = "missing SMILES"
			log_rows.append(log_row)
			continue

		descriptor_values = calculate_descriptors(
			str(smiles_string).strip()
		)

		if descriptor_values is None:
			log_row = dict(base_log)
			log_row["rdkit_status"] = "failed"
			log_row["rdkit_note"] = "RDKit could not parse SMILES"
			log_rows.append(log_row)
			continue

		result_row = {
			"compound_id": compound_id,
			"CAS.No": cas_no,
			"Compound": compound,
			"SMILES": smiles_string,
			"representation": representation_name,
			"structure_source": row["structure_source"],
			"structure_note": row["structure_note"],
			"structure_status": row["structure_status"],
			"structure_type": row["structure_type"]
		}

		result_row.update(descriptor_values)
		results.append(result_row)

		log_row = dict(base_log)
		log_row["rdkit_status"] = "success"
		log_row["rdkit_note"] = "descriptor calculation completed"
		log_rows.append(log_row)

	result_df = pd.DataFrame(results)
	log_df = pd.DataFrame(log_rows)

	return result_df, log_df


############################################################
# Main
############################################################

def main():
	if not os.path.exists(input_path):
		raise FileNotFoundError(f"Input file not found: {input_path}")

	df = pd.read_csv(input_path)

	required_cols = [
		"compound_id",
		"CAS.No",
		"Compound",
		"SMILES_as_reported",
		"SMILES_parent",
		"structure_source",
		"structure_note",
		"structure_status",
		"structure_type"
	]

	missing_cols = [
		col for col in required_cols
		if col not in df.columns
	]

	if missing_cols:
		raise ValueError(
			"Missing required columns: " + ", ".join(missing_cols)
		)

	as_reported_df, as_reported_log = calculate_descriptor_table(
		df=df,
		smiles_col="SMILES_as_reported",
		representation_name="as_reported"
	)

	parent_df, parent_log = calculate_descriptor_table(
		df=df,
		smiles_col="SMILES_parent",
		representation_name="parent"
	)

	log_df = pd.concat(
		[
			as_reported_log,
			parent_log
		],
		ignore_index=True
	)

	os.makedirs(
		os.path.dirname(output_as_reported_path),
		exist_ok=True
	)

	os.makedirs(
		os.path.dirname(output_parent_path),
		exist_ok=True
	)

	os.makedirs(
		os.path.dirname(log_path),
		exist_ok=True
	)

	as_reported_df.to_csv(
		output_as_reported_path,
		index=False
	)

	parent_df.to_csv(
		output_parent_path,
		index=False
	)

	log_df.to_csv(
		log_path,
		index=False
	)

	print("RDKit descriptor calculation complete.")
	print(f"Input compounds: {len(df)}")
	print(f"As-reported descriptors: {len(as_reported_df)}")
	print(f"Parent descriptors: {len(parent_df)}")
	print(f"As-reported output: {output_as_reported_path}")
	print(f"Parent output: {output_parent_path}")
	print(f"Log output: {log_path}")


if __name__ == "__main__":
	main()