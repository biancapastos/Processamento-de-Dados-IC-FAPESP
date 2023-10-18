# Usage
# snakemake

configfile: "config.yaml"

import pandas as pd
from pandas.errors import EmptyDataError
import yaml
import os

numreads_to_integer = config["software"]["numreads_to_integer"]
combats = config["software"]["combats"]
intra_normalizations = config["software"]["intra_normalizations"]
spikes_for_ruvg = config["software"]["spikes"]
ruvg = config["software"]["ruvg"]
inter_normalizations = config["software"]["inter_normalizations"]
qnt_norm = config["software"]["qnt_norm"]
remove_batch_effect = config["software"]["remove_batch_effect"]

estimated_numreads_matrix = "/home/biapastos/snakemake-pipeline/labbces_racs_ic_biancapastos_2023/inputs/Sviridis_v0.2_Matrix_numreads.txt"
batch_file = "/home/biapastos/snakemake-pipeline/labbces_racs_ic_biancapastos_2023/inputs/Batch_table.txt"
gene_lengths = "/home/biapastos/snakemake-pipeline/labbces_racs_ic_biancapastos_2023/inputs/gene_length.txt"
spikes_ruv = "/home/biapastos/snakemake-pipeline/labbces_racs_ic_biancapastos_2023/inputs/ruvg_spikes.txt"
groups_ruv = "/home/biapastos/snakemake-pipeline/labbces_racs_ic_biancapastos_2023/inputs/ruvg_groups.txt"

genotype = "Sviridis_v0.2"

numreads_methods = ['integer', 'combatseq', 'ruvg'] #ruvg tambem
intra_methods = ['cpm', 'tpm', 'fpkm']
inter = ['tmm', 'ctf', 'uq', 'cuf']

integer_or_ruv = ['integer', 'ruvg']
qnt_none = ['qnt', 'none']

rule all:
        input:
            expand("{genotype}/matrices/{genotype}_{numreads_method}_{intra_method}_qnt.txt" , genotype=genotype, numreads_method=numreads_methods, intra_method=intra_methods),
            expand("{genotype}/matrices/{genotype}_{numreads_method}_none_{inter}.txt", genotype=genotype, numreads_method=numreads_methods, inter=inter),
            expand("{genotype}/matrices/rmvbatch_{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}.txt", genotype=genotype, integer_or_ruv=integer_or_ruv, intra_method=intra_methods, qnt_none=qnt_none),
            expand("{genotype}/matrices/removeffect_{genotype}_{integer_or_ruv}_none_{inter}.txt", genotype=genotype, integer_or_ruv=integer_or_ruv, inter=inter),
            expand("{genotype}/matrices/cmbt_{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}.txt", genotype=genotype, integer_or_ruv=integer_or_ruv, intra_method=intra_methods, qnt_none=qnt_none),
            expand("{genotype}/matrices/combat_{genotype}_{integer_or_ruv}_none_{inter}.txt", genotype=genotype, integer_or_ruv=integer_or_ruv, inter=inter)

# MASTER RULE
rule round_to_integer:
        priority: 1
        input:
                expand("{estimated_numreads_matrix}", estimated_numreads_matrix=estimated_numreads_matrix)
        output:
                integer_matrix = "{genotype}/matrices/{genotype}_integer.txt"
        threads: 1
        resources:
                load=1
        log:
                "{genotype}/logs/matrices/{genotype}_integer.log"
        shell:
                """
                python3 {numreads_to_integer} --salmon_numreads {input} --integer_matrix {output.integer_matrix}
                """

rule combatseq_norm:
        input:
                integer_matrix = "{genotype}/matrices/{genotype}_integer.txt",
                batch_table = expand({batch_file}, batch_file=batch_file)
        output:
                combatseq_matrix = "{genotype}/matrices/{genotype}_combatseq.txt"
        threads: 1
        resources:
                load=1
        params:
                method="combatseq"
        log:
                "{genotype}/matrices/{genotype}_combatseq.log"
        shell:
                """
                python3 {combats} --matrix {input.integer_matrix} --batch_table {input.batch_table} --method combatseq --output_matrix {output.combatseq_matrix} 
                """

rule ruvg_norm:
        input:
                integer_matrix = "{genotype}/matrices/{genotype}_integer.txt",
                spikes_file = expand({spikes_ruv}, spikes_ruv=spikes_ruv),
                groups_file = expand({groups_ruv}, groups_ruv=groups_ruv)
        output:
                "{genotype}/matrices/{genotype}_ruvg.txt"
        threads: 1
        resources:
                load=1
        log:
                "{genotype}/{genotype}_ruvg.log"
        shell:
                """
                {ruvg} --integer_matrix {input.integer_matrix} --spikes {input.spikes_file} --groups {input.groups_file} --output_matrix {output}
                """

rule intra_normalizations:
        input:
                matrix = "{genotype}/matrices/{genotype}_{numreads_method}.txt",
                #matrix = expand("{{genotype}}/matrices/{{genotype}}_{numreads_method}.txt", numreads_method=numreads_methods),
                gene_table = expand("{gene_length}", gene_length=gene_lengths)
        output:
                "{genotype}/matrices/{genotype}_{numreads_method}_{intra_method}_none.txt"
        threads: 1
        resources:
                load=1
        params:
                intra_method = "{intra_method}"
        log:
                "{genotype}/logs/matrices/{genotype}_{numreads_method}_{intra_method}.log"
        shell:
                """
                python3 {intra_normalizations} --matrix {input.matrix} --gene_length_table {input.gene_table} --intra_method {params.intra_method} --output_matrix {output}
                """

rule qnt:
        input:
                "{genotype}/matrices/{genotype}_{numreads_method}_{intra_method}_none.txt"
                #expand("{{genotype}}/matrices/{{genotype}}_{numreads_method}_{intra_method}_none.txt", numreads_method=numreads_methods, intra_method=intra_methods)
        output:
                "{genotype}/matrices/{genotype}_{numreads_method}_{intra_method}_qnt.txt"
        threads: 1
        resources:
                load=1
        log:
                "{genotype}/logs/matrices/{genotype}_{numreads_method}_{intra_method}_qnt.log"
        shell:
                """
                {qnt_norm} --matrix {input} --output_matrix {output}
                """

rule inter_normalizations:
        input:
                "{genotype}/matrices/{genotype}_{numreads_method}.txt"
                #expand("{{genotype}}/matrices/{{genotype}}_{numreads_method}.txt", numreads_method=numreads_methods)
        output:
                "{genotype}/matrices/{genotype}_{numreads_method}_none_{inter}.txt"
        threads: 1
        resources:
                load=1
        params:
                norm_method = "{inter}"
        log:
                "{genotype}/logs/matrices/{genotype}_{numreads_method}_none_{inter}.log"
        shell:
                """
                {inter_normalizations} --matrix {input} --norm_method {params.norm_method} --output_matrix {output}
                """

rule removebatch_preprocess_intra:
        input:
                intra_matrices = "{genotype}/matrices/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}.txt",
                #intra_matrices = expand("{{genotype}}/matrices/{{genotype}}_{integer_or_ruv}_{intra_method}_{qnt_none}.txt", integer_or_ruv=integer_or_ruv, intra_method=intra_methods, qnt_none=qnt_none),
                batch_table = expand({batch_file}, batch_file=batch_file)
        output:
                study_matrix = "{genotype}/matrices/studies_batches/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}_std_matx.txt",
                batch_index = "{genotype}/matrices/studies_batches/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}_batchindex.txt"
        threads: 1
        resources:
                load=1
        params:
                method="removebatcheffect"
        #log:
        #        "{genotype}/matrices/{genotype}_{numreads_method}_{intra_none}_{inter_none}_studymatrix.log",
        #        "{genotype}/matrices/{genotype}_{numreads_method}_{intra_none}_{inter_none}_batchindex.log"
        shell:
                """
                python3 {combats} --matrix {input.intra_matrices} --batch_table {input.batch_table} --method {params.method} --output_matrix {output.study_matrix} --batch_index {output.batch_index}
                """

rule removebatch_preprocess_inter:
        input:
                inter_matrices = "{{genotype}}/matrices/{{genotype}}_{integer_or_ruv}_none_{inter}.txt",
                #inter_matrices = expand("{{genotype}}/matrices/{{genotype}}_{integer_or_ruv}_none_{inter}.txt", integer_or_ruv=integer_or_ruv, inter=inter),
                batch_table = expand({batch_file}, batch_file=batch_file)
        output:
                study_matrix = "{genotype}/matrices/studies_batches/studymatrix_{genotype}_{integer_or_ruv}_none_{inter}.txt",
                batch_index = "{genotype}/matrices/studies_batches/batchindex_{genotype}_{integer_or_ruv}_none_{inter}.txt"
        threads: 1
        resources:
                load=1
        params:
                method="removebatcheffect"
        #log:
        #        "{genotype}/matrices/{genotype}_{numreads_method}_{intra_none}_{inter_none}_studymatrix.log",
        #        "{genotype}/matrices/{genotype}_{numreads_method}_{intra_none}_{inter_none}_batchindex.log"
        shell:
                """
                python3 {combats} --matrix {input.inter_matrices} --batch_table {input.batch_table} --method {params.method} --output_matrix {output.study_matrix} --batch_index {output.batch_index}
                """

rule removebatcheffect_intra:
        input:
                study_matrix = "{genotype}/matrices/studies_batches/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}_std_matx.txt",
                batch_index = "{genotype}/matrices/studies_batches/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}_batchindex.txt"
        output:
                removebatch_matrix = "{genotype}/matrices/rmvbatch_{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}.txt"
        threads: 1
        resources:
                load=1
        #log:
        #        "{genotype}/matrices/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}_removebatch.log"
        shell:
                """
                {remove_batch_effect} --matrix {input.study_matrix} --batch_effect {input.batch_index} --output_matrix {output.removebatch_matrix}
                """

rule removebatcheffect_inter:
        input:
                study_matrix = "{genotype}/matrices/studies_batches/studymatrix_{genotype}_{integer_or_ruv}_none_{inter}.txt",
                batch_index = "{genotype}/matrices/studies_batches/batchindex_{genotype}_{integer_or_ruv}_none_{inter}.txt"
        output:
                removebatch_matrix = "{genotype}/matrices/removeffect_{genotype}_{integer_or_ruv}_none_{inter}.txt"
        threads: 1
        resources:
                load=1
        #log:
        #        "{genotype}/matrices/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}_removebatch.log"
        shell:
                """
                {remove_batch_effect} --matrix {input.study_matrix} --batch_effect {input.batch_index} --output_matrix {output.removebatch_matrix}
                """

rule combat_norm_intra:
        input:
                matrix = "{genotype}/matrices/{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}.txt",
                batch_table = expand({batch_file}, batch_file=batch_file)
        output:
                combat_matrix = "{genotype}/matrices/cmbt_{genotype}_{integer_or_ruv}_{intra_method}_{qnt_none}.txt"
        threads: 1
        resources:
                load=1
        params:
                method="combat"
        #log:
         #       "{genotype}/matrices/pycombat/{genotype}_{integer_or_ruv}_{intra_method}_{inter_qnt}_combat.txt"
        shell:
                """
                python3 {combats} --matrix {input.matrix} --batch_table {input.batch_table} --method {params.method} --output_matrix {output.combat_matrix}
                """
rule combat_norm_inter:
        input:
                matrix = "{genotype}/matrices/{genotype}_{integer_or_ruv}_none_{inter}.txt",
                batch_table = expand({batch_file}, batch_file=batch_file)
        output:
                combat_matrix = "{genotype}/matrices/combat_{genotype}_{integer_or_ruv}_none_{inter}.txt"
        threads: 1
        resources:
                load=1
        params:
                method="combat"
        #log:
         #       "{genotype}/matrices/pycombat/{genotype}_{integer_or_ruv}_{intra_method}_{inter_qnt}_combat.txt"
        shell:
                """
                python3 {combats} --matrix {input.matrix} --batch_table {input.batch_table} --method {params.method} --output_matrix {output.combat_matrix}
                """