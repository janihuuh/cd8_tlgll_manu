## Init
import os
import numpy as np
import numpy.random as random
import pandas as pd

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.inference.autotune import auto_tune_scvi_model

from umap import UMAP

import torch
import scanpy as sc
import louvain

import logging
import pickle
from hyperopt import hp


# %matplotlib inline

use_cuda = True
n_epochs_all = None
save_path = ''
show_plot = True
os.chdir("/scratch/cs/csb/projects/lgll/")


## Download samples
FM2311 = CsvDataset(filename='results/scvi/input_files/FM2311.csv', save_path='', sep=',', new_n_genes=False)
fm2477 = CsvDataset(filename='results/scvi/input_files/fm2477.csv', save_path='', sep=',', new_n_genes=False)
FM2481 = CsvDataset(filename='results/scvi/input_files/FM2481.csv', save_path='', sep=',', new_n_genes=False)
fm1189 = CsvDataset(filename='results/scvi/input_files/fm1189.csv', save_path='', sep=',', new_n_genes=False)
fm1890 = CsvDataset(filename='results/scvi/input_files/fm1890.csv', save_path='', sep=',', new_n_genes=False)
fm1338 = CsvDataset(filename='results/scvi/input_files/fm1338.csv', save_path='', sep=',', new_n_genes=False)
FM1656_17 = CsvDataset(filename='results/scvi/input_files/FM1656_17.csv', save_path='', sep=',', new_n_genes=False)
FM1656_18 = CsvDataset(filename='results/scvi/input_files/FM1656_18.csv', save_path='', sep=',', new_n_genes=False)
FM1382 = CsvDataset(filename='results/scvi/input_files/FM1382.csv', save_path='', sep=',', new_n_genes=False)
t805_15 = CsvDataset(filename='results/scvi/input_files/805_15.csv', save_path='', sep=',', new_n_genes=False)
t805_18 = CsvDataset(filename='results/scvi/input_files/805_18.csv', save_path='', sep=',', new_n_genes=False)


all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs = [FM2311.X,
                                               fm2477.X,
                                               FM2481.X,
                                               fm1189.X,
                                               fm1890.X,
                                               fm1338.X,
                                               FM1656_17.X,
                                               FM1656_18.X,
                                               FM1382.X,
                                               t805_15.X,
                                               t805_18.X])


## Train, save and fin
vae      = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels, n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')
trainer  = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=100)
torch.save(trainer.model.state_dict(), 'results/scvi/results/lgll_oneshot.pkl')

## Sample posterior to get latent representation and save those embeddings
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))

latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

np.savetxt("results/scvi/results/lgll_oneshot_latent.csv", latent, delimiter=",")
np.savetxt("results/scvi/results/lgll_oneshot_indices.csv", batch_indices, delimiter=",")
