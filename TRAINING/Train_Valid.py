from keras.models import Sequential
from keras.layers import Dense, Activation
import matplotlib.pyplot as plt
import ast
import numpy as np
import keras

file_res = open("TRAINING\\results.txt", "r")
file_par = open("TRAINING\\training_params.txt", "r")
lines_res = file_res.read()
lines_par = file_par.read()
results = list(filter(('').__ne__, lines_res.split('\n')))
params = list(filter(('').__ne__, lines_par.split('\n')))
ones = []
train = []
k = 0
for i in range(0, len(results)):
    if results.index(results[i]) % 3 == 0:
        ones = ones + list(ast.literal_eval(results[i]))
        train = train + list(ast.literal_eval(params[int(i/3)]))
        if len(ones) < len(train):
            ones = ones + [1]*(len(train)-len(ones))
print(len(ones), len(train))     
trainX = np.array(train)
trainY = np.array(ones)


# weight_1 = np.array(ast.literal_eval(arr.split('array(')[1].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# bias_1 = np.array(ast.literal_eval(arr.split('array(')[2].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# weight_2 = np.array(ast.literal_eval(arr.split('array(')[3].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# bias_2 = np.array(ast.literal_eval(arr.split('array(')[4].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# weight_3 = np.array(ast.literal_eval(arr.split('array(')[5].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)
# bias_3 = np.array(ast.literal_eval(arr.split('array(')[6].split('dtype=')[0].replace("],\n", "]").replace("\n", "")), dtype=np.float32)

model = Sequential()
model.add(Dense(32, input_dim=15)) 
model.add(Activation('relu'))
model.add(Dense(16)) 
model.add(Activation('relu'))
model.add(Dense(1))
model.add(Activation('sigmoid'))


# j = 1
# for layer in [model.layers[0], model.layers[2], model.layers[4]]:
#     layer.set_weights([eval(f"weight_{j}")[0], eval(f"bias_{j}")[0]])
#     j = j + 1

checkpoint_path = "c:\\Users\\Pedro\\DCYFR\\weights_checkpoint.keras"
checkpoint = keras.callbacks.ModelCheckpoint(filepath=checkpoint_path, save_freq="epoch", monitor='accuracy', mode='max', save_best_only=True)
model.compile(loss=keras.losses.BinaryFocalCrossentropy(gamma=2.0, alpha=0.25), optimizer=keras.optimizers.Adamax(learning_rate=0.001), metrics=['accuracy']) #
history = model.fit(trainX, trainY, epochs=100, batch_size=80, verbose=1, validation_split=0.25, callbacks=[checkpoint])

keras.models.load_model(checkpoint_path)

weighty_file = open("weights_3.txt", "w")
weighty_file.write(str(model.get_weights()))
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()