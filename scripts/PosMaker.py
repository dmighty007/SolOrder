import numpy as np
from tqdm import tqdm

class PosArray:
    def __init__(self, biased, step = 1):
        self.step = step
        self.biased = biased
        self.n_frames = biased.trajectory.n_frames
        self.n_atom = biased.select_atoms("backbone").n_atoms
        self.pos_array = np.zeros(((self.n_frames//step), self.n_atom, 3))
        self.extract_positions()
        self.shift_and_scale_positions()
        self.flatten_data()

    def extract_positions(self):
        for i, frame in tqdm(enumerate(self.biased.trajectory[::self.step ]), total=self.n_frames // self.step ):
            self.pos_array[i] = self.biased.select_atoms("backbone").positions

    def shift_and_scale_positions(self):
        shifted_pos = np.array(list(map(lambda x: x-x[0], self.pos_array)))
        #scale = abs(shifted_pos).mean()
        scale = 1
        self.scale_shifed_pos = shifted_pos / scale

    def flatten_data(self):
        self.flat_data = np.array(list(map(lambda x : x.ravel(), self.scale_shifed_pos)))