HUMANCHROMOSOMES=["10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "1", "20", "21", "22", "2", "3", "4", "5", "6", "7", "8", "9", "X"]
MOUSECHROMOSOMES=["10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "1", "2", "3", "4", "5", "6", "7", "8", "9", "X"]

rule cutMaps:
	input:
		"data/HiC/{dataset}.nonnormalized/{chromosome}.txt"
	output:
		"data/{dataset}.nonnormalized/{chromosome}.txt.mod"
	shell:
		"for file in {input}; do tail -n +3 $file | cut -f2- > {output}; done"

rule normalizeMaps:
	input:
		hmaps=expand("data/40kb.nonnormalized/chr{chr}.nij.combined.40kb.matrix.LA.NULLS.txt.mod", chr=HUMANCHROMOSOMES),
		mmaps=expand("data/mESC.nonnormalized/nij.chr{chr}.LA.txt.mod", chr=MOUSECHROMOSOMES)
	output:
		"data/{dataset}.nonnormalized/{chromosome}.txt.mod.dmat"
	shell:
		"scripts/Normalizer.py {input.hmaps} {input.mmaps}"

rule skipHeader:
	input:
		"data/genefiles/{organism}/{organism}HomMartExp.txt"
	output:
		"data/genefiles/{organism}HomMartExpWhoTit.txt"
	shell:
		"tail -n +2 {input} > {output}"

rule makeAnnotationFile:
	input:
		"data/genefiles/{organism}HomMartExpWhoTit.txt"
	output:
		"data/genefiles/{organism}Ensembl.txt"
	shell:
		"scripts/AnFileMaker.py {input} {output}"

rule makeHomologyTable:
	input:
		"data/genefiles/humanHomMartExpWhoTit.txt",
		"data/genefiles/mouseHomMartExpWhoTit.txt"
	output:
		"data/genefiles/HomologyHS_MM_Ensembl.csv"
	shell:
		"scripts/MakeHomList.py {input} {output};"
		"uniq {output} > tmp.csv;"
		"mv tmp.csv {output}"

rule buildHumanGraphs:
	input:
		dmat=expand("data/HiC/40kb.nonnormalized/chr{chr}.nij.combined.40kb.matrix.LA.NULLS.txt.mod.dmat", chr=HUMANCHROMOSOMES),
		txt="data/genefiles/humanEnsembl.txt",
		csv="data/genefiles/HomologyHS_MM_Ensembl.csv"
	output:
		"GraphhumanSelfNormIntraEnsSeqOrd.ml",
		"GraphhumanSelfNormIntraEns3D.ml"
	shell:
		"scripts/SeqGraphMaker.py -s 40000 {input.txt} {input.csv} GraphhumanSelfNormIntraEnsSeqOrd.ml {input.dmat} -o 10 11 12 13 14 15 16 17 18 19 1 20 21 22 2 3 4 5 6 7 8 9 X;"
		"scripts/ParseToGraphml.py -s 40000 {input.txt} {input.csv} GraphhumanSelfNormIntraEns3D.ml {input.dmat} -o 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 1 1 20 20 21 21 22 22 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 X X"

rule buildMouseGraphs:
	input:
		dmat=expand("data/HiC/mESC.nonnormalized/nij.chr{chr}.LA.txt.mod.dmat", chr=MOUSECHROMOSOMES),
		txt="data/genefiles/mouseEnsembl.txt",
		csv="data/genefiles/HomologyHS_MM_Ensembl.csv"
	output:
		"GraphmouseSelfNormIntraEnsSeqOrd.ml",
		"GraphmouseSelfNormIntraEns3D.ml"
	shell:
		"scripts/SeqGraphMaker.py -s 40000 {input.txt} {input.csv} GraphmouseSelfNormIntraEnsSeqOrd.ml {input.dmat} -o 10 11 12 13 14 15 16 17 18 19 1 2 3 4 5 6 7 8 9 X;"
		"scripts/ParseToGraphml.py -s 40000 {input.txt} {input.csv} GraphmouseSelfNormIntraEns3D.ml {input.dmat} -o 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 X X"

rule findTeams:
	input:
		graph=expand("Graph{org}SelfNormIntraEns{graphtype}.ml", org=["human", "mouse"], graphtype=["SeqOrd", "3D"])
	output:
		"selfnorm_ensemble_delta_d{delta}_{graphtype}.csv"
	params:
		delt="{delta}"
	shell:
		"scripts/graph_teams.py -d {params.delt} {input} > {output}"
