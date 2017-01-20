# atac-seq-pipeline
ATAC-seq analysis pipeline

## Dependencies
- Python
- [bowtie][1]
- [samtools][2]
- [bedtools][3]
- [MACS2][4]
- [numpy][5]
- [IDR][6]

## align-reads.py


## call-atac.py

```
Usage: python run-idr.py [options] <expt prefix> <output dir> <genome sz>

Options:
  -m   : Use parameters manually specified in script.
         This option will ignore following arguments.
```

## run-idr.py

```
Usage: python call-atac.py [options] <input dir> <expt prefix> <output dir> <genome sz>

Options:
  -m   : Use parameters manually specified in script.
         This option will ignore following arguments.
```



[1]:	https://github.com/BenLangmead/bowtie
[2]:	http://www.htslib.org/doc/samtools.html
[3]:	http://bedtools.readthedocs.io/en/latest/
[4]:	https://github.com/taoliu/MACS
[5]:	http://www.numpy.org/
[6]:	https://github.com/nboley/idr
