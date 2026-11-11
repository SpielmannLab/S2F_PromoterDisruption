#!/usr/bin/env python
import os
from pathlib import Path

import pandas as pd
from alphagenome_research.model import dna_model
from dotenv import load_dotenv
from huggingface_hub import login

print("Now loaded all the libraries")

# Log in to huggingface (Set the MY_TOKEN_PATH env variable in the shell)
load_dotenv(dotenv_path=Path(os.environ["MY_TOKEN_PATH"]))
login()
print("Now logged in to HuggingFace")

# Set up a dna_model (which is a server)
model = dna_model.create_from_huggingface("all_folds")
print("Started the dna_model")

# Gather all supported ontologies (biological samples)
cage_ontologies = model.output_metadata().cage
rnaseq_ontologies = model.output_metadata().rna_seq

rnaseq_ontologies.to_csv(
    "./input/alphagenome_available_ontologies_rnaseq.tsv", sep="\t"
)
cage_ontologies.to_csv("./input/alphagenome_available_ontologies_cage.tsv", sep="\t")

# Now merge them to predict only common ones
shared_ontologies = pd.merge(
    rnaseq_ontologies,
    cage_ontologies,
    on=("strand", "ontology_curie", "biosample_name", "biosample_type"),
    how="inner",
    suffixes=("_rnaseq", "_cage"),
)

shared_ontologies.to_csv(
    "./input/alphagenome_available_ontologies_shared.tsv", sep="\t"
)
