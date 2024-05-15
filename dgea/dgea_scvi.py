import numpy as np
import pandas as pd
import random
import anndata as ad

import torch
import scvi

model_to_load = torch.load("model.pt")

