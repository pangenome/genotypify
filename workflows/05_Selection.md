# Selection

## Variables

```shell
dir_base=/lizardfs/guarracino/genotypify
```

## bmws

Run it:

```shell
mkdir -p $dir_base/amylase/bmws && $dir_base/amylase/bmws
conda activate /lizardfs/guarracino/condatools/bmws/0.2.1

bmws analyze $dir_base/data/amy.vcf.gz $dir_base/data/amy.meta -d diploid -l 4.5 -g 30 -n 10000 -t -o $dir_base/amylase/bmws/amy.out
```

Plot the result:

```R
library(tidyverse)
library(cowplot)

s_mean <- read_lines("amy.out") %>%
  str_split(pattern = "\t") %>%
  .[[1]] %>%
  .[7] %>%
  parse_double()

s_trajectory <- read_lines("amy.out") %>%
  str_split(pattern = "\t") %>%
  .[[1]] %>%
  .[-(1:9)] %>% 
  parse_double() %>%
  tibble(s=.) %>%
  transmute(generation=1:387, time=(387-generation)*30, s)

# Plot trajectory per generation
s_trajectory %>%
  ggplot(aes(x=generation, y=s)) +
  geom_line()

# Plot trajectiory per time
s_trajectory %>%
  ggplot(aes(x=-time/1000, y=s)) +
  geom_line() +
  geom_hline(yintercept = s_mean, linetype=3) +
  scale_x_continuous(breaks = -4:0*3, labels = 4:0*3) +
  scale_y_continuous(breaks=0:3 * 0.02) +
  ylim(c(0, 0.06)) +
  xlab("kyr BP") +
  theme_cowplot()
```

Parameter estimation:

```shell
python3 $dir_base/scripts/amy1.py $dir_base/data/amy.txt $dir_base/amylase/bmws 1000

# Estimate s
mkdir -p $dir_base/amylase/bmws/bootstrap
conda activate /lizardfs/guarracino/condatools/bmws/0.2.1
for SEED in {0..999}; do
    echo "Running simulation for SEED=${SEED}..."
    python $dir_base/scripts/amy2.py $dir_base/amylase/bmws $dir_base/amylase/bmws/bootstrap ${SEED} > $dir_base/amylase/bmws/bootstrap/amy.s_hat.${SEED}.log 2>&1
done
```

Plot the result:

```R
library(tidyverse)
library(cowplot)
library(ggplot)
library(ggtext)

s_hat <- read_tsv("/home/guarracino/Desktop/Garrison/amylase/bmws/amy.s_hat.txt", col_names = "s") %>%
  mutate(generation=388:1, time=(388-generation)*30)
s_hat_mean=mean(s_hat$s)

paths <- read_tsv("/home/guarracino/Desktop/Garrison/amylase/bmws/amy.paths.tsv", col_names = FALSE) %>%                      
  pivot_longer(cols = 1:389, names_to = "name", values_to = "p") %>%
  mutate(generation=390-parse_number(name), time=(389-generation)*30) %>%
  group_by(time, generation) %>%
  summarise(p_mean=mean(p), p_lower = quantile(p, 0.025), p_higher=quantile(p, 0.975))

s_ci <- read_tsv(str_c("/home/guarracino/Desktop/Garrison/amylase/bmws/bootstrap/amy.s_hat.", 0:999, ".txt"), col_names = "s") %>%
  mutate(generation=rep(388:1, 1000), time=(387-generation)*30) %>%
  group_by(time, generation) %>%
  summarise(s_mean=mean(s), s_lower = quantile(s, 0.025), s_higher=quantile(s, 0.975))

s_plot <- s_hat %>%
  ggplot(aes(x=-time/1000)) +
  geom_line(aes(y=s), color="blue", size=1) +
  #geom_line(data=s_ci, mapping = aes(y=s_mean)) +
  geom_ribbon(data=s_ci, aes(ymin=s_lower, ymax=s_higher), fill="lightblue", alpha=0.5) +
  annotate(geom = "text", x=-2.7, y=s_hat_mean*1.3, label="bar(s)", parse=TRUE, size=5) +
  annotate(geom = "text", x=-1, y=s_hat_mean*1.3, label=str_c("=", round(s_hat_mean, 3)), size=5) +
  geom_hline(yintercept = s_hat_mean, linetype=3) +
  scale_x_continuous(breaks = -4:0*3, labels = 4:0*3) +
  scale_y_continuous(breaks=0:3 * 0.02) +
  xlab("kyr BP") +
  ylab("selection<br>coefficient") +
  theme_cowplot() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color="black", size=1),
        axis.title.y = ggtext::element_markdown(lineheight = 1.2))
s_plot

p_plot <- paths %>%
  ggplot(aes(x=-time/500)) +
  geom_line(aes(y=p_mean), color="blue", size=1) +
  geom_ribbon(aes(ymin=p_lower, ymax=p_higher), fill="lightblue", alpha=0.5) +
  scale_x_continuous(breaks = -4:0*3, labels = 4:0*3) +
  xlab("kyr BP") +
  #ylab(expression(str_c("frequency of\n", italic("dup"), " haplotypes"))) +
  ylab("frequency of<br>*dup* haplotypes") +
  theme_cowplot() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color="black", size=1),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(lineheight = 1.2))
combined_plot <- plot_grid(p_plot, s_plot , nrow = 2, rel_heights  = c(1, 1.12), align = "v")
combined_plot
```


#### Prepare the envirorment and parameter estimation

/home/bioinfo26/Amylase_project/bmws/amy1.py


#### Estimate s 

```shell
for SEED in {0..999}; do
    echo "Running simulation for SEED=${SEED}..."
    source ~/.bashrc
    conda activate bmws
    python /home/bioinfo26/Amylase_project/bmws/amy.2.py ${SEED} > /home/bioinfo26/Amylase_project/bmws/bootstrap/amy.s_hat.${SEED}.log 2>&1
done
```

