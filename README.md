# Processamento-de-Dados-IC-FAPESP
Processamento de arquivos de matrizes, normalização dos dados com Python e R e geração de arquivos de redes. Automatização dos códigos com Snakemake.

## Normalização dos dados
- Métodos de Normalização Intra-amostrais
  - CPM, TPM e FPKM
- Métodos de Normalização Inter-amostrais
  - TMM, CTF, UQ, CUF e QNT
- Métodos de Remoção de Efeitos de Lote
  - ComBat, ComBat-Seq e RemoveBatchEffect
- Método de Remoção de Variações Indesejadas
  - RUV, mais precisamente aplicado o método RUVg que utiliza spikes (genes de controle negativo)

### Arquivos
- /inputs
  - Contém arquivos de entrada utilizados nos códigos como a matriz NumReads gerada pelo Salmon, arquivo com o comprimento do gene a ser utilizado nas normalizações TPM e FPKM, arquivo com os batches a serem utilizados para os métodos de remoção de efeito de lote e arquivos com os spikes e grupos utilizados no RUVg.
- /scripts
  - Contém os códigos em Python ou R que utilizam bibliotecas/pacotes para realizar a normalização dos dados contidos nas matrizes de entrada e geram novos arquivos com uma matriz de valores normalizados.
- Snakefile.smk
  - Código que implementa o snakemake com as regras que fazem a automatização de todos os arquivos de matrizes a serem gerados.
