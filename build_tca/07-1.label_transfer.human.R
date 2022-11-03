options(stringsAsFactors = FALSE)
library(here)
library(magrittr)
library(tidyverse)
library(patchwork)
library(mlr3verse)
source(here("scripts/R/utils.R"))
source(here("scripts/R/00.theme.R"))
source(here("scripts/R/MeasureMultiRocAuc.R"))


OUTDIR <- here("results/07.label_transfer.human")
safe_mkdir(OUTDIR)

#### 1. build dataset ####
dataset <- readRDS(here("results/06.cNMF_human/hTCA_mc.all.cnmf.hvg2000/K80/hTCA.AUCell.rds"))

files <- list.files(here("results/03.manual_annotation/"), full.names = T)
files <- files[grepl("*human*", files)]

celltype.df <- lapply(files, function(xx) {
  sel.cols <- c("Cell_type_manual", "Cell_type_predicted")
  read.table(xx, header = T, sep = "\t", row.names = 1)[, sel.cols]
}) %>% do.call(rbind, .)

dataset$Cell_type_manual <- mapvalues.2(rownames(dataset), from = rownames(celltype.df), to = celltype.df$Cell_type_manual)
dataset$Cell_type_manual <- ifelse(dataset$Cell_type_manual == "Unknown", NA, dataset$Cell_type_manual)
dataset.labeled <- subset(dataset, !is.na(Cell_type_manual))
dataset.labeled$Cell_type_manual %<>% factor()
dataset.unlabeled <- subset(dataset, is.na(Cell_type_manual))
dataset.unlabeled$Cell_type_manual <- NULL

saveRDS(dataset.labeled, file.path(OUTDIR, "hTCA.AUCell.labeled.rds"))
saveRDS(dataset.unlabeled, file.path(OUTDIR, "hTCA.AUCell.unlabeled.rds"))

#### 2. model benchmark ####
dataset.labeled <- readRDS(file.path(OUTDIR, "hTCA.AUCell.labeled.rds"))

task = as_task_classif(dataset.labeled, target = "Cell_type_manual", id = "hTCA")
design = benchmark_grid(
  tasks = task,
  learners = lrns(c("classif.xgboost", "classif.randomForest", "classif.svm"),
                  predict_type = "prob", predict_sets = c("test")),
  resamplings = rsmps("cv", folds = 5)
)
bmr = benchmark(design)
measures = list(
  msr("classif.acc", predict_sets = "test", id = "acc_test"),
  msr("mclassif.mauc_aunu", predict_sets = "test", id = "auc_test")
)
scores = bmr$score(measures) %>% as.data.frame()

## plots
p1 <- scores %>%
  mutate(learner_id = sub("classif.", "", learner_id, fixed = T)) %>%
  ggplot(aes(learner_id, acc_test, color=learner_id)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  labs(y="Accuracy", x="Model", title = "Model performance") +
  theme(legend.position = "none")

p2 <- scores %>%
  mutate(learner_id = sub("classif.", "", learner_id, fixed = T)) %>%
  ggplot(aes(learner_id, auc_test, color=learner_id)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  labs(y="AUC (ROC)", x="Model", title = "Model performance") +
  theme(legend.position = "none")

p1 + p2
ggsave(file.path(OUTDIR, "model.benchmark.pdf"), width = 8, height = 4)

#### 3. label transfer ####

dataset.labeled <- readRDS(file.path(OUTDIR, "hTCA.AUCell.labeled.rds"))
dataset.unlabeled <- readRDS(file.path(OUTDIR, "hTCA.AUCell.unlabeled.rds"))
dataset <- rbind(dataset.labeled[-ncol(dataset.labeled)], dataset.unlabeled)
task = as_task_classif(dataset.labeled, target = "Cell_type_manual", id = "hTCA")

## SVM
learner = lrn("classif.svm")
learner$predict_type = "prob"
learner$train(task)

pred.label <- learner$predict_newdata(newdata = dataset)
dataset$pred.label <- pred.label$response
dataset$purity <- apply(pred.label$prob, 1, max)

saveRDS(learner, file.path(OUTDIR, "hTCA.manual.SVM.rds"))
saveRDS(dataset[, c("pred.label", "purity")], file.path(OUTDIR, "hTCA.predicted_labels.rds"))


