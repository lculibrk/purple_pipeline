import yaml
import glob
import re

chrs = [
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    "X"
]

with open("pog_libs.yaml") as f:
    pog_libs = yaml.load(f, Loader = yaml.FullLoader)

#print(pog_libs)

pog_libs = {k:v for (k, v) in pog_libs.items() if re.match("POG", k)}

reverso = {}
for case in pog_libs.keys():
    reverso[case] = {}
    #print(pog_libs[case])
    for sample in pog_libs[case].keys():
#        print(sample)
         glob_res = glob.glob(pog_libs[case][sample]["merge"] + "/*dupsFlagged.bam")
         if not glob_res:
            continue
         reverso[case][pog_libs[case][sample]["name"]] = glob_res[0]

## Build input list

in_list = []

id_dict = pog_libs

cases=id_dict.keys()

blacklist=["A79702","A79703","A79704"]

#print(cases)

comps={}
bams={}
for case in cases:
    #print(id_dict[case])
    if id_dict[case].get("Normal_1") and id_dict[case].get("Diseased_1"):
#        print(case)
        norm = id_dict[case]["Normal_1"]
        bams[norm["name"]] = glob.glob(norm["merge"] + "/*dupsFlagged.bam")
        if len(bams[norm["name"]]) == 0:
            continue
        lib_iter = 1
        bams[norm["name"]] = bams[norm["name"]][0]
        tums = []
#        print(id_dict[case])
#        print("Diseased_" + str(lib_iter))
#        print(id_dict[case].get("Diseased_" + str(lib_iter)))
        while id_dict[case].get("Diseased_" + str(lib_iter)):
#            print(id_dict[case]["Diseased_" + str(lib_iter)])
            bam_glob = glob.glob(id_dict[case]["Diseased_" + str(lib_iter)]["merge"] + "/*dupsFlagged.bam")
            if len(bam_glob) == 0:
                lib_iter += 1
                continue
            if os.path.islink(bam_glob[0]) and not os.path.exists(bam_glob[0]):
                lib_iter += 1
                continue
            if id_dict[case]["Diseased_" + str(lib_iter)]["name"] in blacklist:
                lib_iter += 1
                continue
            bams[id_dict[case]["Diseased_" + str(lib_iter)]["name"]] = bam_glob[0]
            tums.append(id_dict[case]["Diseased_" + str(lib_iter)]["name"] + "_" + norm["name"])
            lib_iter += 1
        comps[case] = tums
        #print(tums)
        #libs = list(id_dict[case].keys())
        #libs = re.
        #print(libs)
in_list = []
tools = [
    #"cobalt",
    "amber"#,
    #"purple"]
]
amber_list = []
cobalt_list = []
purple_list = []
print(comps)
for comp in comps.keys():
    for compstring in comps[comp]:
        compstring = compstring.split("_")
        tumor = compstring[0]
        normal = compstring[1]
        amber_list.extend(expand("amber/{case}/{tumor}_{normal}/", case = comp, tumor = tumor, normal = normal))
        cobalt_list.extend(expand("cobalt/{case}/{tumor}_{normal}/", case = comp, tumor = tumor, normal = normal))
        purple_list.extend(expand("purple/{case}/{tumor}_{normal}/", case = comp, tumor = tumor, normal = normal))
#        in_list.extend(expand("{tool}/{case}/{tumor}_{normal}/{normal}.amber.snp.vcf.gz", tool = tools, case = comp, tumor = tumor, normal = normal))



counter = 0
chunk = 1
amber_dict = {}
chunk_list = []
for item in amber_list:
## Parse item
    tumor = re.sub(".*/([^_]*)_([^_]*)", "\\1", item)
    normal = re.sub(".*/([^_]*)_([^_]*)/", "\\2", item)
    s_item = re.sub("amber/(.*)", "\\1", item)
    if counter == 99 or item == amber_list[-1]:
        counter = 0
        amber_dict["chunk" + str(chunk)] = chunk_list
        chunk = chunk + 1
        chunk_list = []
    chunk_list.append("chunks/amber/chunk" + str(chunk) + "/" + s_item + tumor + ".amber.baf.vcf.gz")
    chunk_list.append("chunks/amber/chunk" + str(chunk) + "/" + s_item + normal + ".amber.snp.vcf.gz")
    counter = counter + 1
print(amber_dict.keys())
amber_list = amber_dict.values()

counter = 0
chunk = 1
cobalt_dict = {}
chunk_list = []
for item in cobalt_list:
## Parse item
    tumor = re.sub(".*/([^_]*)_([^_]*)", "\\1", item)
    normal = re.sub(".*/([^_]*)_([^_]*)/", "\\2", item)
    s_item = re.sub("cobalt/(.*)", "\\1", item)
    if counter == 99 or item == cobalt_list[-1]:
        counter = 0
        cobalt_dict["chunk" + str(chunk)] = chunk_list
        chunk = chunk + 1
        chunk_list = []
    chunk_list.append("chunks/cobalt/chunk" + str(chunk) + "/" + s_item + tumor + ".cobalt.ratio.pcf")
    chunk_list.append("chunks/cobalt/chunk" + str(chunk) + "/" + s_item + normal + ".cobalt.ratio.pcf")
    counter = counter + 1
print(cobalt_dict.keys())
cobalt_list = cobalt_dict.values()

counter = 0
chunk = 1
purple_dict = {}
chunk_list = []
for item in purple_list:
## Parse item
    tumor = re.sub(".*/([^_]*)_([^_]*)", "\\1", item)
    normal = re.sub(".*/([^_]*)_([^_]*)/", "\\2", item)
    s_item = re.sub("purple/(.*)", "\\1", item)
    if counter == 99 or item == purple_list[-1]:
        counter = 0
        purple_dict["chunk" + str(chunk)] = chunk_list
        chunk = chunk + 1
        chunk_list = []
    chunk_list.append("chunks/purple/chunk" + str(chunk) + "/" + s_item)
    counter = counter + 1

purple_list = list(purple_dict.values())
purple_list = [item for sublist in purple_list for item in sublist]
chunks = cobalt_dict.keys()
chunks = [int(re.sub("chunk([0-9]+)", "\\1", ch)) for ch in chunks]
chunk = "chunks/chunk" + str(max(chunks))
chunk_output = chunk + "/chunk_complete.txt"


#purple_list = [f for f in purple_list if re.match("chunk1$", f)]

wildcard_constraints:
    chunk = "chunk[0-9]*",
    tumor = "[A-Z][0-9]*",
    normal = "[A-Z][0-9]*",
    case = "[^\/]*"


#print(in_list)

def amber_chunk_function(chunk):
    chunk = int(re.sub("chunk([0-9]+)", "\\1", chunk))
    chunk = "chunks/amber/chunk" + str(chunk - 1)
    return(chunk + "/chunk_complete.txt")

def cobalt_chunk_function(chunk):
    chunk = int(re.sub("chunk([0-9]+)", "\\1", chunk))
    chunk = "chunks/cobalt/chunk" + str(chunk - 1)
    return(chunk + "/chunk_complete.txt")

localrules: all, init_chunk, run_chunk_amber, run_chunk_cobalt

rule all:
    input:
#        amber_list,
#        cobalt_list,
        purple_list

rule init_chunk:
    output:
        "chunks/amber/chunk0/chunk_complete.txt",
        "chunks/cobalt/chunk0/chunk_complete.txt"
    resources: cpus=1, mem_mb=7900
    shell:
        "touch {output}"

rule run_chunk_amber:
    input:
        lambda wildcards: amber_chunk_function(wildcards.chunk),
        lambda wildcards: amber_dict[wildcards.chunk]
    output:
        "chunks/amber/{chunk}/chunk_complete.txt"
    resources: cpus=1, mem_mb=7900
    shell:
        "touch {output}"

rule run_chunk_cobalt:
    input:
        lambda wildcards: cobalt_chunk_function(wildcards.chunk),
        lambda wildcards: cobalt_dict[wildcards.chunk]
    output:
        "chunks/cobalt/{chunk}/chunk_complete.txt"
    resources: cpus=1, mem_mb=7900
    shell:
        "touch {output}"

rule amber:
    input:
        lambda wildcards: bams[wildcards.normal],
        lambda wildcards: bams[wildcards.tumor],
        lambda wildcards: amber_chunk_function(wildcards.chunk)
    output:
        "chunks/amber/{chunk}/{case}/{tumor}_{normal}/{tumor}.amber.baf.vcf.gz",
        "chunks/amber/{chunk}/{case}/{tumor}_{normal}/{normal}.amber.snp.vcf.gz"
    resources: cpus=5, mem_mb=40000
    benchmark:
        "benchmarks/amber/{chunk}/{case}/{tumor}_{normal}/benchmark.txt"
    shell:
        """
        java -jar software/amber-3.5.jar \
            -reference {wildcards.normal} \
            -reference_bam {input[0]} \
            -tumor {wildcards.tumor} \
            -tumor_bam {input[1]} \
            -output_dir chunks/amber/{wildcards.chunk}/{wildcards.case}/{wildcards.tumor}_{wildcards.normal}/ \
            -threads {resources.cpus} \
            -loci GermlineHetPon.hg19.vcf.gz \
            -validation_stringency LENIENT
        """

rule cobalt:
    input:
        lambda wildcards: bams[wildcards.normal],
        lambda wildcards: bams[wildcards.tumor],
        lambda wildcards: cobalt_chunk_function(wildcards.chunk)
    output:
        "chunks/cobalt/{chunk}/{case}/{tumor}_{normal}/{tumor}.cobalt.ratio.pcf",
        "chunks/cobalt/{chunk}/{case}/{tumor}_{normal}/{normal}.cobalt.ratio.pcf"
    resources: cpus=8, mem_mb=40000
    benchmark:
        "benchmarks/cobalt/{chunk}/{case}/{tumor}_{normal}/benchmark.txt"
    shell:
        """
        java -jar software/cobalt-1.11.jar \
            -reference {wildcards.normal} \
            -reference_bam {input[0]} \
            -tumor {wildcards.tumor} \
            -tumor_bam {input[1]} \
            -output_dir chunks/cobalt/{wildcards.chunk}/{wildcards.case}/{wildcards.tumor}_{wildcards.normal}/ \
            -threads {resources.cpus} \
            -gc_profile cobalt_GC_profile.hg19.1000bp.cnp \
            -validation_stringency LENIENT
        """
rule purple:
    input:
        "chunks/cobalt/{chunk}/{case}/{tumor}_{normal}/{tumor}.cobalt.ratio.pcf",
        "chunks/cobalt/{chunk}/{case}/{tumor}_{normal}/{normal}.cobalt.ratio.pcf",
        "chunks/amber/{chunk}/{case}/{tumor}_{normal}/{tumor}.amber.baf.vcf.gz",
        "chunks/amber/{chunk}/{case}/{tumor}_{normal}/{normal}.amber.snp.vcf.gz"
    output:
        directory("chunks/purple/{chunk}/{case}/{tumor}_{normal}")
    resources: cpus=8, mem_mb=40000
    benchmark:
        "benchmarks/purple/{chunk}/{case}/{tumor}_{normal}/benchmark.txt"
    shell:
        """
        java -jar software/purple-2.51.jar \
            -reference {wildcards.normal} \
            -tumor {wildcards.tumor} \
            -output_dir {output} \
            -amber chunks/amber/{wildcards.chunk}/{wildcards.case}/{wildcards.tumor}_{wildcards.normal}/ \
            -cobalt chunks/cobalt/{wildcards.chunk}/{wildcards.case}/{wildcards.tumor}_{wildcards.normal}/ \
            -gc_profile cobalt_GC_profile.hg19.1000bp.cnp \
            -ref_genome resources/hg19a.fasta
        """
