# Oral microbiome of patients under suspicion of gastric alterations: a prospective, bicentric cross-sectional study

To reproduce, please build the Docker image:

```bash
docker build -t microbiome-and-gastric-alterations:0.0.1 .
```

Then run each script:

```bash
docker run --rm -i -v ${PWD}:/home/rstudio/ microbiome-and-gastric-alterations:0.0.1 Rscript analysis.R
```

and

```bash
docker run --rm -i -v ${PWD}:/home/rstudio/ microbiome-and-gastric-alterations:0.0.1 Rscript analysis-stratified.R
```