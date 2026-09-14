import os
import re
from pathlib import Path

configfile: "config.yaml"
R1_PATTERNS = [
    r'^(.*)(_1)(\.fq|\.fastq)(\.gz)?$',
    r'^(.*)(_R1_001)(\.fq|\.fastq)(\.gz)?$',
    r'^(.*)(_R1)(\.fq|\.fastq)(\.gz)?$',
    r'^(.*)(\.1)(\.fq|\.fastq)(\.gz)?$',
]

R2_PATTERNS = [
    r'^(.*)(_2)(\.fq|\.fastq)(\.gz)?$',
    r'^(.*)(_R2_001)(\.fq|\.fastq)(\.gz)?$',
    r'^(.*)(_R2)(\.fq|\.fastq)(\.gz)?$',
    r'^(.*)(\.2)(\.fq|\.fastq)(\.gz)?$',
]

PE_REGEXES = [re.compile(p) for p in (*R1_PATTERNS, *R2_PATTERNS)]

FASTQ_EXTENSIONS = (
    ".fq",
    ".fastq",
    ".fq.gz",
    ".fastq.gz",
)


def get_prefix(filename):
    """Return sample prefix for PE or SE filename."""
    for rx in PE_REGEXES:
        m = rx.match(filename)
        if m:
            return m.group(1)

    for ext in sorted(FASTQ_EXTENSIONS, key=len, reverse=True):
        if filename.endswith(ext):
            return filename[:-len(ext)]

    return None


def discover_samples(sampledir_path):
    sampledir = Path(sampledir_path)
    samples = {}

    for root, _, files in os.walk(sampledir):
        rel_dir = Path(root).relative_to(sampledir)

        for filename in files:
            prefix = get_prefix(filename)
            if prefix is None:
                continue

            if rel_dir.name == prefix:
                sample_id = str(rel_dir)
            else:
                sample_id = (
                    prefix
                    if rel_dir == Path(".")
                    else str(rel_dir / prefix)
                )

            samples.setdefault(sample_id, []).append(
                str(Path(root) / filename)
            )

    for files in samples.values():
        files.sort()

    return samples


# Initialize parameters
SAMPLEDIR = Path(config["sampledir"])
SAMPLES_DICT = discover_samples(SAMPLEDIR)

PE_SAMPLES = [s for s, files in SAMPLES_DICT.items() if len(files) == 2]
SE_SAMPLES = [s for s, files in SAMPLES_DICT.items() if len(files) == 1]
ALL_SAMPLES = list(SAMPLES_DICT.keys())

print(f"Discovered {len(ALL_SAMPLES)} samples total.")
print(f"  Paired-End: {len(PE_SAMPLES)}")
print(f"  Single-End: {len(SE_SAMPLES)}")

MINLEN = config.get("trimmer", {}).get("minlen", 120)
OUTDIR = config["outdir"]
QCTYPE = config.get("qctype", "trimmed")
if QCTYPE == "": 
   QCTYPE = "trimmed"

rule all:
   input:
        # Raw FastQC
        fastqc_pe_1=expand(OUTDIR + "/fastqc/{ID}_1_fastqc.zip", ID=PE_SAMPLES),
        fastqc_pe_2=expand(OUTDIR + "/fastqc/{ID}_2_fastqc.zip", ID=PE_SAMPLES),
        fastqc_se=expand(OUTDIR + "/fastqc/{ID}_1_fastqc.zip", ID=SE_SAMPLES),
        
        # Trimmed FastQC
        fastqcp_pe_1=expand(OUTDIR + "/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip", ID=PE_SAMPLES),
        fastqcp_pe_2=expand(OUTDIR + "/fastqc/{ID}/2."+QCTYPE+"_fastqc.zip", ID=PE_SAMPLES),
        fastqcp_se=expand(OUTDIR + "/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip", ID=SE_SAMPLES),
        
        # Centralized Links
        OGlink_pe_1=expand(OUTDIR + "/fastqc/OG/{ID}/1_fastqc.zip", ID=PE_SAMPLES),
        OGlink_pe_2=expand(OUTDIR + "/fastqc/OG/{ID}/2_fastqc.zip", ID=PE_SAMPLES),
        OGlink_se=expand(OUTDIR + "/fastqc/OG/{ID}/1_fastqc.zip", ID=SE_SAMPLES),
        
        QCTYPElink_pe_1=expand(OUTDIR + "/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip", ID=PE_SAMPLES),
        QCTYPElink_pe_2=expand(OUTDIR + "/fastqc/"+QCTYPE+"/{ID}/2_fastqc.zip", ID=PE_SAMPLES),
        QCTYPElink_se=expand(OUTDIR + "/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip", ID=SE_SAMPLES),
   output:
        multiqc1=OUTDIR+"/fastqc/OG/QC_report_multiqc_report.html",
        multiqc2=OUTDIR+"/fastqc/"+QCTYPE +"/QC_report_multiqc_report.html"
   shell:
       """
       multiqc {OUTDIR}/fastqc/OG --interactive --title "QC_report" -f -o {OUTDIR}/fastqc/OG
       multiqc {OUTDIR}/fastqc/{QCTYPE} --interactive --title "QC_report" -f -o {OUTDIR}/fastqc/{QCTYPE}
       echo Done
       """

# --- PAIRED-END RULES ---

rule fastqc_pe:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
        f2=lambda wildcards: SAMPLES_DICT[wildcards.ID][1],
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/{ID}_2_fastqc.zip"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
        
        tmp2=$(mktemp -d)
        fastqc {input.f2} -o $tmp2
        mv $tmp2/*_fastqc.zip {output.fastqc2}
        rm -rf $tmp2
       """

rule fastqcp_pe:
   input:
        f1=OUTDIR +"/" +QCTYPE+"/{ID}_1.fq.gz",
        f2=OUTDIR +"/" +QCTYPE+"/{ID}_2.fq.gz"
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/{ID}/2."+QCTYPE+"_fastqc.zip"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
        
        tmp2=$(mktemp -d)
        fastqc {input.f2} -o $tmp2
        mv $tmp2/*_fastqc.zip {output.fastqc2}
        rm -rf $tmp2
       """

rule fastp_pe:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
        f2=lambda wildcards: SAMPLES_DICT[wildcards.ID][1],
   output:
        fastqc1=OUTDIR+"/"+QCTYPE+"/{ID}_1.fq",
        fastqc2=OUTDIR+"/"+QCTYPE+"/{ID}_2.fq",
   params:
        adapters=config.get("adapters", "adapters.fa"),
   log:
        stdout=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stdout.txt",
        stderr=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stderr.txt",
   threads: 8
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
       mkdir -p $(dirname {output.fastqc1})/log

       if [ "{QCTYPE}" == "fastp" ];
       then
               if [ "{params.adapters}" != "" ] && [ -f "{params.adapters}" ];
               then
                       ADAPTER="--adapter_fasta {params.adapters}" 
               else
                       ADAPTER=""
               fi
               
               fastp --overrepresentation_analysis $ADAPTER --detect_adapter_for_pe --trim_poly_g \
                  --cut_mean_quality 20 --trim_poly_x --length_required 30 --thread {threads} \
                  --cut_window_size 4 --cut_front --cut_tail  \
                  --html $(dirname {output.fastqc1})/{QCTYPE}.html  \
                  --json $(dirname {output.fastqc1})/{QCTYPE}.json \
                  -i {input.f1}  -I {input.f2} -o {output.fastqc1}.tmp -O {output.fastqc2}.tmp > {log.stdout}  2>{log.stderr}

               cutadapt -a 'A{{10}}' -a 'T{{10}}' -A 'A{{10}}' -A 'T{{10}}' \
                  -o $(dirname {output.fastqc1})/fastp1c \
                  -p $(dirname {output.fastqc1})/fastp2c \
                  {output.fastqc1}.tmp {output.fastqc2}.tmp >> {log.stdout} 2>>{log.stderr}
                  
               rm {output.fastqc1}.tmp {output.fastqc2}.tmp
               
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i $(dirname {output.fastqc1})/fastp1c -o {output.fastqc1} -v >> {log.stdout} 2>>{log.stderr}
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i $(dirname {output.fastqc1})/fastp2c -o {output.fastqc2} -v >> {log.stdout} 2>>{log.stderr}
               
               rm $(dirname {output.fastqc1})/fastp1c $(dirname {output.fastqc1})/fastp2c
       else
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f1}) -o {output.fastqc1} -v > {log.stdout} 2>{log.stderr}
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f2}) -o {output.fastqc2} -v >> {log.stdout} 2>>{log.stderr}
       fi
       """

rule gzip_pe:
   input:
        fastqc1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq",
        fastqc2=OUTDIR + "/"+QCTYPE+"/{ID}_2.fq"
   output:
        gz1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq.gz",
        gz2=OUTDIR + "/"+QCTYPE+"/{ID}_2.fq.gz"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
       """
       gzip {input.fastqc1}
       gzip {input.fastqc2}
       """

rule link_pe:
   input:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip",
        fastqc2=OUTDIR+"/fastqc/{ID}_2_fastqc.zip",
        f1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
        f2=OUTDIR+"/fastqc/{ID}/2."+QCTYPE+"_fastqc.zip"
   output:
        OGln1=OUTDIR+"/fastqc/OG/{ID}/1_fastqc.zip",
        OGln2=OUTDIR+"/fastqc/OG/{ID}/2_fastqc.zip",
        FPln1=OUTDIR+"/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip",
        FPln2=OUTDIR+"/fastqc/"+QCTYPE+"/{ID}/2_fastqc.zip",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in PE_SAMPLES) if PE_SAMPLES else "$^"
   shell:
        """
        mkdir -p $(dirname {output.OGln1})
        mkdir -p $(dirname {output.FPln1})
        ln -sf $(readlink -f {input.fastqc1}) {output.OGln1}
        ln -sf $(readlink -f {input.fastqc2}) {output.OGln2}
        ln -sf $(readlink -f {input.f1}) {output.FPln1}
        ln -sf $(readlink -f {input.f2}) {output.FPln2}
        """

# --- SINGLE-END RULES ---

rule fastqc_se:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip"
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
       """

rule fastqcp_se:
   input:
        f1=OUTDIR +"/" +QCTYPE+"/{ID}_1.fq.gz",
   output:
        fastqc1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
        mkdir -p $(dirname {output.fastqc1})
        
        tmp1=$(mktemp -d)
        fastqc {input.f1} -o $tmp1
        mv $tmp1/*_fastqc.zip {output.fastqc1}
        rm -rf $tmp1
       """

rule fastp_se:
   input:
        f1=lambda wildcards: SAMPLES_DICT[wildcards.ID][0],
   output:
        fastqc1=OUTDIR+"/"+QCTYPE+"/{ID}_1.fq",
   params:
        adapters=config.get("adapters", "adapters.fa"),
   log:
        stdout=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stdout.txt",
        stderr=OUTDIR +"/"+QCTYPE+ "/{ID}/log/stderr.txt",
   threads: 8
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
       mkdir -p $(dirname {output.fastqc1})/log

       if [ "{QCTYPE}" == "fastp" ];
       then
               if [ "{params.adapters}" != "" ] && [ -f "{params.adapters}" ];
               then
                       ADAPTER="--adapter_fasta {params.adapters}" 
               else
                       ADAPTER=""
               fi
               
               fastp --overrepresentation_analysis $ADAPTER --trim_poly_g \
                  --cut_mean_quality 20 --trim_poly_x --length_required 30 --thread {threads} \
                  --cut_window_size 4 --cut_front --cut_tail  \
                  --html $(dirname {output.fastqc1})/{QCTYPE}.html  \
                  --json $(dirname {output.fastqc1})/{QCTYPE}.json \
                  -i {input.f1} -o {output.fastqc1}.tmp > {log.stdout}  2>{log.stderr}

               cutadapt -a 'A{{10}}' -a 'T{{10}}' \
                  -o $(dirname {output.fastqc1})/fastp1c \
                  {output.fastqc1}.tmp >> {log.stdout} 2>>{log.stderr}
                  
               rm {output.fastqc1}.tmp
               
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i $(dirname {output.fastqc1})/fastp1c -o {output.fastqc1} -v >> {log.stdout} 2>>{log.stderr}
               
               rm $(dirname {output.fastqc1})/fastp1c
       else
               /mnt/lustre/RDS-live/downing/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l {MINLEN} -i <(zcat {input.f1}) -o {output.fastqc1} -v > {log.stdout} 2>{log.stderr}
       fi
       """

rule gzip_se:
   input:
        fastqc1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq",
   output:
        gz1=OUTDIR + "/"+QCTYPE+"/{ID}_1.fq.gz",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
       """
       gzip {input.fastqc1}
       """

rule link_se:
   input:
        fastqc1=OUTDIR+"/fastqc/{ID}_1_fastqc.zip",
        f1=OUTDIR+"/fastqc/{ID}/1."+QCTYPE+"_fastqc.zip",
   output:
        OGln1=OUTDIR+"/fastqc/OG/{ID}/1_fastqc.zip",
        FPln1=OUTDIR+"/fastqc/"+QCTYPE+"/{ID}/1_fastqc.zip",
   wildcard_constraints:
        ID="|".join(re.escape(s) for s in SE_SAMPLES) if SE_SAMPLES else "$^"
   shell:
        """
        mkdir -p $(dirname {output.OGln1})
        mkdir -p $(dirname {output.FPln1})
        ln -sf $(readlink -f {input.fastqc1}) {output.OGln1}
        ln -sf $(readlink -f {input.f1}) {output.FPln1}
        """
