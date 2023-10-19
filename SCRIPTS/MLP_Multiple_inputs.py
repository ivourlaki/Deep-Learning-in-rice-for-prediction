#!/usr/bin/env python
# coding: utf-8


# Import modules
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt
from sklearn.preprocessing import scale
import sys


# Import tensorflow
import tensorflow as tf
print('tensorflow: %s' % tf.__version__)
from tensorflow import keras
print('keras: %s' % keras.__version__)
import keras_tuner as kt
print('kerastuner: %s' % kt.__version__)



############################## 1. Insert phenotype and markers #################################################

############### arguments import by the terminal script in which we run the python script #####################

y_test = pd.read_csv('/final_partitions_pheno/pheno_{0}_partitions_{1}_test.csv'.format(trait,d),
                         header=None)
y_train = pd.read_csv('/final_partitions_pheno/pheno_{0}_partitions_{1}_train.csv'.format(trait,d),
                          header=None)


X_train1 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_train_final_kernel_snps.csv'.format(trait,d),
                           header=None)
X_test1 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_test_final_kernel_snps.csv'.format(trait,d),
                          header=None)

X_train2 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_train_final_kernel_del.csv'.format(trait,d),
                           header=None)
X_test2 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_test_final_kernel_del.csv'.format(trait,d),
                          header=None)

X_train3 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_train_final_kernel_dup.csv'.format(trait,d),
                           header=None)
X_test3 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_test_final_kernel_dup.csv'.format(trait,d),
                          header=None)

X_train4 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_train_final_kernel_inv.csv'.format(trait,d),
                           header=None)
X_test4 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_test_final_kernel_inv.csv'.format(trait,d),
                          header=None)

X_train5 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_train_final_kernel_mitedtx.csv'.format(trait,d),
                           header=None)
X_test5 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_test_final_kernel_mitedtx.csv'.format(trait,d),
                          header=None)

X_train6 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_train_final_kernel_rlxrix.csv'.format(trait,d),
                           header=None)
X_test6 = pd.read_csv('/final_partitions_pca/trait_{0}_partitions_{1}_test_final_kernel_rlxrix.csv'.format(trait,d),
                          header=None)
  

print(X_train1.shape,X_train2.shape,X_train3.shape,X_train4.shape,X_train5.shape,X_train6.shape,
     X_test1.shape,X_test2.shape,X_test3.shape,X_test4.shape,X_test5.shape,X_test6.shape)


############################ 2. Data preparation: ###########################################################

########################### for this case data scaling has been implemented in R ###########################
################### Test here which type of trait you use #######################################################

# check how many unique values you have
y = y_train
test = y.unique()

print(test)


if  len(test) == 2 : 
    NUM_CLASSES = "binary" # or just 0
    print("you have to recode in 0 and 1")
    
    # 1. check your two values 
    min(test)
    max(test)
    
   
    # 2. replace each one by 0, 1. The minimum by 0 and the maximum by 1. 
    y_train1=y_train.replace(min(test) , 0)
    y_train2=y_train1.replace(max(test) , 1)
    y_train=y_train2
    
    y_test1=y_test.replace(min(test) , 0)
    y_test2=y_test1.replace(max(test) , 1)
    y_test=y_test2
    
    # 3. check the new dataframes
    print(y_train, y_test)

    # 4. check if you have actually only 0 and 1 values. 
    y1 = y_train ; y2 = y_test
    test1 = y1.unique(); test2 = y2.unique()
    print(test1, test2)
    print(NUM_CLASSES)

    

elif isinstance(test[1], float) :
    print(isinstance(test, float))
    NUM_CLASSES = 0
    print("Quantitative trait")
    print(NUM_CLASSES)

else: 
    
    # we dont use this case
    NUM_CLASSES= max(test)
    print("Multiclass trait")
    NUM_CLASSES= NUM_CLASSES + 1
    print(NUM_CLASSES)

y_train=np.squeeze(y_train)
y_test=np.squeeze(y_test)


# they are all the same, I separate them though.

nSNP1 = X_train1.shape[1] 
nSNP2 = X_train1.shape[1] 
nSNP3 = X_train1.shape[1] 
nSNP4 = X_train1.shape[1] 
nSNP5 = X_train1.shape[1] 
nSNP6 = X_train1.shape[1] 
print(nSNP1,nSNP2,nSNP3,nSNP4,nSNP5,nSNP6)



############################ 3. Functional Model for MLP multiple inputs #######################################


if 'model' in locals(): del model # deletes current model if exists

from tensorflow.keras import layers
from tensorflow.keras import Model
from keras_tuner import HyperModel


class build_model(HyperModel):

    def __init__(self, input_dim, num_classes):
            # input and output sizes, if needed, are defined
            self.input_dim = input_dim
            self.num_classes = num_classes


    def build(self,hp):

        inputA = layers.Input(shape=(self.input_dim))
        inputB = layers.Input(shape=(self.input_dim))
        inputC = layers.Input(shape=(self.input_dim))
        inputD = layers.Input(shape=(self.input_dim))
        inputE = layers.Input(shape=(self.input_dim))
        inputF = layers.Input(shape=(self.input_dim))

        
        unit=hp.Choice("units", values=[16,38,64,128])
        kernels = hp.Choice('l', values=[0.001,0.01,0.1])
        x1 = layers.Dense(units=unit,activation='relu', 
                               kernel_regularizer=keras.regularizers.L1L2(l1=kernels,l2=kernels),
                                bias_regularizer=keras.regularizers.L2(kernels),
                                activity_regularizer=keras.regularizers.L2(kernels))(inputA)

        x2 = layers.Dense(units=unit,activation='relu', 
                               kernel_regularizer=keras.regularizers.L1L2(l1=kernels,l2=kernels),
                                bias_regularizer=keras.regularizers.L2(kernels),
                                activity_regularizer=keras.regularizers.L2(kernels))(inputB)

        x3 = layers.Dense(units=unit,activation='relu', 
                               kernel_regularizer=keras.regularizers.L1L2(l1=kernels,l2=kernels),
                                bias_regularizer=keras.regularizers.L2(kernels),
                                activity_regularizer=keras.regularizers.L2(kernels))(inputC)

        x4 = layers.Dense(units=unit,activation='relu', 
                               kernel_regularizer=keras.regularizers.L1L2(l1=kernels,l2=kernels),
                                bias_regularizer=keras.regularizers.L2(kernels),
                                activity_regularizer=keras.regularizers.L2(kernels))(inputD)

        x5 = layers.Dense(units=unit,activation='relu', 
                               kernel_regularizer=keras.regularizers.L1L2(l1=kernels,l2=kernels),
                                bias_regularizer=keras.regularizers.L2(kernels),
                                activity_regularizer=keras.regularizers.L2(kernels))(inputE)

        x6 = layers.Dense(units=unit,activation='relu', 
                               kernel_regularizer=keras.regularizers.L1L2(l1=kernels,l2=kernels),
                                bias_regularizer=keras.regularizers.L2(kernels),
                                activity_regularizer=keras.regularizers.L2(kernels))(inputF)

        
        ## Concatenate the layers here
        
        x =  layers.Concatenate(axis=1)([x1,x2,x3,x4,x5,x6])

        activation=hp.Choice("activation", values=["relu", "tanh","linear"])

        for i in range(hp.Int('n_layers', 0,4)):
            kernels= hp.Choice(f'kernels_{i}_units', values=[0.001,0.01,0.1])
            x = layers.Dense(units=hp.Choice(f'dense_{i}_units',values=[2,4,8,16]), 
                         activation=activation, kernel_regularizer=keras.regularizers.L1L2(l1=kernels,l2=kernels),
                                bias_regularizer=keras.regularizers.L2(kernels),
                                activity_regularizer=keras.regularizers.L2(kernels) )(x)

        x = layers.Dropout(hp.Float('rate', min_value=0.0, max_value=0.3, step=0.05))(x)

        if self.num_classes==0:    # regression
            activation= "linear"
            units=1
          
            x = layers.Dense(units, activation=activation)(x)

            model = Model(inputs=[inputA, inputB,inputC,inputD,inputE,inputF], outputs=x)

            model.compile(optimizer=hp.Choice('optimizer', values=['Adam', 'RMSprop','SGD']),
                    loss='mse', 
                    metrics=["mse"])

        elif self.num_classes == "binary": # binary
                activation = "sigmoid"
                units = 1
               

                x = layers.Dense(units, activation=activation)(x)

                model = Model(inputs=[inputA, inputB,inputC,inputD,inputE,inputF], outputs=x)

                model.compile(optimizer=hp.Choice('optimizer', values=['Adam', 'RMSprop','SGD']),
                loss='binary_crossentropy',
                metrics=["binary_accuracy"])

        else:                  # multiple classes
                activation = "softmax"
                units = self.num_classes
            
                x = layers.Dense(units, activation=activation)(x)

                model = Model(inputs=[inputA, inputB,inputC,inputD,inputE,inputF], outputs=x)

                model.compile(optimizer=hp.Choice('optimizer', values=['Adam', 'RMSprop','SGD']),
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'])

        return model

  
    
hypermodel = build_model(input_dim=nSNP1, num_classes=NUM_CLASSES)


############################   4. Choose the tuner ####################################################################
#  You may choose from RandomSearch, BayesianOptimization and Hyperband.


if NUM_CLASSES == 0 :
    tuner = kt.Hyperband(hypermodel,
                     objective="val_mse",
                     max_epochs=15,
                     hyperband_iterations=20,
                     overwrite=True, directory="LOG_DIR")

   # search summary
    tuner.search_space_summary()

else :
    tuner = kt.Hyperband(hypermodel,
                     objective="val_binary_accuracy",
                     max_epochs=15,
                     hyperband_iterations=20,
                     overwrite=True,directory="LOG_DIR")

 # search summary
    tuner.search_space_summary()



#########################   5.start search ################################################################################


if NUM_CLASSES == 0 :

    stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience = 5, 
                                        restore_best_weights = True)
    tensor=tf.keras.callbacks.TensorBoard(log_dir="LOG_DIR1")
    tuner.search([X_train1,X_train2,X_train3,X_train4,X_train5,X_train6],  y_train, epochs=25, validation_split=0.2,
            # Use the TensorBoard callback.
            # The logs will be write to "/tmp/tb_logs".
                callbacks=[stop_early,tensor] 
    )


else :
# monitor it could be accuracy
    stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss',patience = 5, 
                                        restore_best_weights = True)
    
    tensor=tf.keras.callbacks.TensorBoard(log_dir="LOG_DIR1")
    tuner.search([X_train1,X_train2,X_train3,X_train4,X_train5,X_train6], y_train, epochs=25, validation_split=0.2,
            # Use the TensorBoard callback.
            # The logs will be write to "/tmp/tb_logs".
                callbacks=[stop_early,tensor])



#########################   6.Show the results! ############################################################################

tuner.results_summary()


###########################  7.Retrieve the best model or hyperparmeters ####################

best_weight = tuner.get_best_models(num_models=1)[0]




best_hps=tuner.get_best_hyperparameters(num_trials=1)[0]
print(best_hps.values)
model = tuner.hypermodel.build(best_hps) 
model.summary()




######################## 8.save the best param ###############################################################################
bestparam=best_hps.values

file = open(f"partition_{d}_bestparam.txt","w")
 
for key, value in bestparam.items(): 
 
 file.write('%s:%s\n' % (key, value))
 
file.close()



###################### 9. Retrain the model ##################################################################################

stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5 ,
                                            restore_best_weights = True)

batch_size=32
epoch=25
best_model = tuner.hypermodel.build(best_hps)
df1=pd.DataFrame()

X_train1, X_val, y_train1, y_val = train_test_split(X_train, y_train, test_size=0.20)

for x in range(1):    
    test_model= best_model
    test_model.fit([X_train1,X_train2,X_train3,X_train4,X_train5,X_train6], y_train, epochs=epoch,validation_split=0.2,
                            batch_size=batch_size, callbacks=[stop_early]) 
    
    model_new = test_model
   
   

    if NUM_CLASSES == 0:



        y_hat = model_new.predict([X_test1,X_test2,X_test3,X_test4,X_test5,X_test6])
        test_scores = model_new.evaluate([X_test1,X_test2,X_test3,X_test4,X_test5,X_test6], y_test, verbose=2)
        print("Test loss:", test_scores[0])


         # correlation btw predicted and observed
        mse = tf.keras.metrics.mean_squared_error(
        y_test,np.squeeze(y_hat))
        coef = np.corrcoef(y_test,y_hat[:,0])[0,1]
        print('\nCorr obs vs pred =',mse)
        corr2= pd.DataFrame([mse,coef], index=["loss","coef"])
        df=pd.DataFrame(corr2)
        df1=pd.concat([df1,df], axis=1, ignore_index=True)



    elif NUM_CLASSES == "binary":
        #a.  predict the classes


        y_hat = model_new.predict([X_test1,X_test2,X_test3,X_test4,X_test5,X_test6])
        test_scores = model_new.evaluate([X_test1,X_test2,X_test3,X_test4,X_test5,X_test6], y_test, verbose=2)
        print("Test loss:", test_scores[0])
        
        y_hat01_1= np.where(y_hat > 0.5, 1,0)
        print(y_hat01_1.shape)
        coef = np.corrcoef(y_test,y_hat01_1[:,0])[0,1]
        bce=tf.keras.metrics.BinaryCrossentropy()
        bce1 =bce(y_test, np.squeeze(y_hat)).numpy()
        print(bce1)
        acc =tf.keras.metrics.BinaryAccuracy()
        acc1=acc(y_test, np.squeeze(y_hat)).numpy()
        print(acc1)
        corr2= pd.DataFrame([acc1, bce1, coef], index=["accuracy","loss","coef"])
        df=pd.DataFrame(corr2)
        df1=pd.concat([df1,df], axis=1, ignore_index=True)




####################### 10. Save the evaluation metric values ########################################################################

df1.to_csv(f"eval__metrics_partition_{d}.csv")       





