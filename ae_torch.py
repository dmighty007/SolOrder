import numpy as np
import tqdm.auto as tqdm
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

class AE(nn.Module):
    
    def __init__(self, input_dim,train_data, val_data, hidden1 = 128, hidden2 = 40,
                 code = 2, learning_rate=0.00001,batch_size=64, epochs=50,
                 shuffle=True, thresh = 0.0001):
        
        """
        input_dim : flattened input vector length
        hidden1 : node number of hidden layer 1
        hidden2 : node number of hidden layer 2
        code : dimension of latent space
        learning_rate : name suggests
        thresh : thresh to compare while earlystopping
        train_data : trainning dataset
        val_data : validation dataset
        """
        
        super(AE, self).__init__()
        
        self.input_dim = input_dim
        self.batch_size = batch_size
        self.epochs = epochs
        self.shuffle = shuffle
        self.thresh = thresh
        self.learning_rate = learning_rate
        self.train_data = train_data
        self.val_data = val_data
        
        torch.manual_seed(1)
        np.random.seed(1)
        # Encoder block
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden1),
            nn.ReLU(),
            nn.Linear(hidden1, hidden2),
            nn.ReLU(),
            nn.Linear(hidden2, code),
        )
        # Decoder block
        self.decoder = nn.Sequential(
            nn.Linear(code, hidden2),
            nn.ReLU(),
            nn.Linear(hidden2, hidden1),
            nn.ReLU(),
            nn.Linear(hidden1, input_dim),
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded
    
    def fit(self):
        
        if torch.cuda.is_available():
            self.to('cuda')
        train_dataset = TensorDataset(torch.Tensor(self.train_data))
        train_dataloader = DataLoader(train_dataset, batch_size=self.batch_size, shuffle=self.shuffle)

        val_dataset = TensorDataset(torch.Tensor(self.val_data))
        val_dataloader = DataLoader(val_dataset, batch_size=self.batch_size, shuffle=self.shuffle)
        
        criterion = nn.MSELoss()
        optimizer = optim.Adam(self.parameters(), lr=self.learning_rate)
        best_val_loss = float('inf')
        prv = float('inf')
        stop_counter = 0
        for epoch in range(self.epochs):
            train_loss = 0
            for x_batch, in train_dataloader:
                x_batch = x_batch.to('cuda') if torch.cuda.is_available() else x_batch
                optimizer.zero_grad()
                decoded = self(x_batch)
                loss = criterion(decoded, x_batch)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()

            val_loss = 0
            for x_batch, in val_dataloader:
                x_batch = x_batch.to('cuda') if torch.cuda.is_available() else x_batch
                decoded = self(x_batch)
                loss = criterion(decoded, x_batch)
                val_loss += loss.item()
            val_loss = val_loss / len(val_dataloader)
            if abs(val_loss - prv) < self.thresh:
                break
            elif val_loss < best_val_loss:
                best_val_loss = val_loss
                stop_counter = 0
                # You can save the best model using torch.save(self.state_dict(), 'best_model.pt')
           
            else:
                stop_counter += 1
                if stop_counter >= 2:
                    print('Early stopping')
                    break
                    
            prv = val_loss

            print(f'Epoch: {epoch + 1}/{self.epochs} | '
                  f'Train Loss: {train_loss / len(train_dataloader):.5f} | '
                  f'Val Loss: {val_loss:.5f} | ')
