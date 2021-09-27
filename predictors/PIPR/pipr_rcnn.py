"""
    ---------- Original Work this file is based on: ----------
    Published in Bioinformatics journal featuring ISMB/ECCB 2019
    Title: 'Multifaceted Protein-Protein Interaction Prediction Based on Siamese Residual RCNN'
    Authors: Chen, Muhao and Ju, Chelsea and Zhou, Guangyu and Chen, Xuelu and Zhang, Tianran and Chang, Kai-Wei and Zaniolo, Carlo and Wang, Wei
    Journal: Bioinformatics
    Volume: 35
    Number: 14
    Pages: i305-i314
    Year: 2019
    Month: 07
    Publisher: Oxford University Press
    DOI: http://dx.doi.org/10.1093/bioinformatics/btz328
    git: https://github.com/muhaochen/seq_ppi
    
    ---------- This file ----------
    This pipr_rcnn.py file is a modification from the original git file seq_ppi/binary/model/lasagna/rcnn.py
    Main modifications include a change of command-line argument usage for execution and a choice of cross-validation 
    or a single train/test split. Prediction probabilities of each interaction in test data are also saved to file.
    Author: Eric Arezza
    Last Updated: March 9, 2021
    
    Description:
        RCNN approach to binary classification of protein-protein interaction prediction.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys, os, argparse
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import tensorflow as tf
import keras.backend.tensorflow_backend as KTF
import numpy as np
import pickle

from datetime import datetime
from time import time
from embeddings.seq2tensor import s2t
from tqdm import tqdm

from keras.models import Model
from keras.layers import Dense, Bidirectional, Input, CuDNNGRU, LeakyReLU
from keras.layers.merge import concatenate, multiply
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import MaxPooling1D, GlobalAveragePooling1D
from keras.optimizers import Adam,  RMSprop
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, average_precision_score

if 'embeddings' not in sys.path:
    sys.path.append('../../../embeddings')

# Description of command-line usage
describe_help = 'CUDA_VISIBLE_DEVICES=0 python pipr_rcnn.py sequencesFile.fasta trainFile.tsv testFile.tsv'
parser = argparse.ArgumentParser(description=describe_help)
parser.add_argument('sequences', help='Path to file containing protein sequences (.fasta)', type=str, nargs=1)
parser.add_argument('train', help='Path to file containing binary protein interactions for training (.tsv)', type=str, nargs=1)
parser.add_argument('test', help='Path to file containing binary protein interactions for testing (.tsv)', type=str, nargs=1)
parser.add_argument('-r','--results', help='Path to file to store results', type=str, nargs=1, required=False)
parser.add_argument('-p','--predictions', help='Path to file to store the predictions', type=str, required=False)
parser.add_argument('-l', '--label_index', help='Label index (int)', type=int, nargs=1, required=False)
parser.add_argument('-m', '--mbedding', help='Embedding (int), 0=embeddings/default_onehot.txt, 1=embeddings/string_vec5.txt, 2=embeddings/CTCoding_onehot.txt, 3=embeddings/vec7_CTC.txt', 
                    type=int, nargs=1, required=False, choices=[0,1,2,3])
parser.add_argument('-d', '--dimensions', help='Hidden dimensions (int)', type=int, nargs=1, required=False)
parser.add_argument('-e', '--epochs', help='Epochs (int)', type=int, nargs=1, required=False)
parser.add_argument('-a', '--seq_size', help='Amino acids/sequence length (int)', type=int, nargs=1, required=False)
# parser.add_argument('-save', '--saveModel', help='Save model', action='store_true', default=False)
parser.add_argument('-save', '--saveModel', help='Path where the model is to be saved.', required=False, type=str)
parser.add_argument('-load','--loadModel', help='Path to pre-trained model', default='', type=str, required=False)
parser.add_argument('-c', '--cpu', dest='cpu', help='Use only CPU', action='store_true', default=False)
#parser.set_defaults(results=os.getcwd()+'/results/'+datetime.now().strftime("%d-%m-%Y/")+datetime.now().strftime("%H-%M-%S-results.txt"),
#                    label_index=2, mbedding=0, dimensions=25, epochs=50)
args = parser.parse_args()
# Set defaults for command-line arguments
if args.label_index is None:
    label_index = 2
else:
    label_index = args.label_index[0]
if args.mbedding is None:
    use_emb = 3
else:
    use_emb = args.mbedding[0]
if args.dimensions is None:
    hidden_dim = 50
else:
    hidden_dim = args.dimensions[0]
if args.epochs is None:
    n_epochs = 100
else:
    n_epochs = args.epochs[0]
if args.seq_size is None:
    SEQ_SIZE = 2000
else:
    SEQ_SIZE = args.seq_size[0]
if args.loadModel == '':
    pretrained = None
else:
    pretrained = args.loadModel

TRAIN_FILE = args.train[0]
TEST_FILE = args.test[0]
CROSS_VALIDATE = False
if TRAIN_FILE == TEST_FILE:
    CROSS_VALIDATE = True

rst_file = args.results

ID2SEQ_FILE = args.sequences[0]

EMB_FILES = ['embeddings/default_onehot.txt', 'embeddings/string_vec5.txt', 'embeddings/CTCoding_onehot.txt', 'embeddings/vec7_CTC.txt']
SEQ2T = s2t(EMB_FILES[use_emb])

print("\n---Using the following---\nSequences File: {}\nTraining File: {}\nTesting File: {}\nResults File: {}".format(ID2SEQ_FILE, TRAIN_FILE, TEST_FILE, rst_file))
print("Label index: {}\nEmbedding: {} - {}\nHidden Dimensions: {}\nEpochs: {}\n".format(label_index, use_emb, EMB_FILES[use_emb], hidden_dim, n_epochs))
print('Save model: {}\nLoad model: {}'.format(args.saveModel, pretrained))
DIM = SEQ2T.dim
#SEQ_SIZE = 2000
CLASS_MAP = {'0':1,'1':0}
print("Class map:", CLASS_MAP)

def get_session(gpu_fraction=0.5):
    '''Assume that you have 6GB of GPU memory and want to allocate ~3GB'''

    num_threads = int(os.environ.get('OMP_NUM_THREADS'))
    gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=gpu_fraction)

    if num_threads:
        return tf.Session(config=tf.ConfigProto(
            gpu_options=gpu_options, intra_op_parallelism_threads=num_threads))
    else:
        return tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))

def get_interaction_data(file, sid1_index, sid2_index, id2index, seqs, seq_array):
    
    sid = 0
    id2_aid = {}
    max_data = -1
    limit_data = max_data > 0
    skip_head = True
    count = 0
    raw_data = []
    
    # Assigns proteins an ID number and add sequences to an array indexed same as protein
    for line in tqdm(open(file)):
        if skip_head:
            skip_head = False
            continue
        line = line.rstrip('\n').strip('\r').split('\t')
        # Ensures all proteins in sequences file have an index number
        if id2index.get(line[sid1_index]) is None or id2index.get(line[sid2_index]) is None:
            continue
        # Otherwise, generate new index number
        if id2_aid.get(line[sid1_index]) is None:
            id2_aid[line[sid1_index]] = sid
            sid += 1
            seq_array.append(seqs[id2index[line[sid1_index]]])
        line[sid1_index] = id2_aid[line[sid1_index]]
        if id2_aid.get(line[sid2_index]) is None:
            id2_aid[line[sid2_index]] = sid
            sid += 1
            seq_array.append(seqs[id2index[line[sid2_index]]])
        line[sid2_index] = id2_aid[line[sid2_index]]
        raw_data.append(line)
        if limit_data:
            count += 1
            if count >= max_data:
                break
    print('Raw data interactions:', len(raw_data))
    
    #len_m_seq = np.array([len(line.split()) for line in seq_array])
    #avg_m_seq = int(np.average(len_m_seq)) + 1
    #max_m_seq = max(len_m_seq)
    #print (avg_m_seq, max_m_seq)

    return raw_data, seq_array, id2_aid

def create_class_labels(raw_data, class_map, label_index):
    class_labels = np.zeros((len(raw_data), 2))
    for i in range(len(raw_data)):
        class_labels[i][class_map[raw_data[i][label_index]]] = 1.
    return class_labels

def build_model():
    x = None
    seq_input1 = Input(shape=(SEQ_SIZE, DIM), name='seq1')
    seq_input2 = Input(shape=(SEQ_SIZE, DIM), name='seq2')
    l1=Conv1D(hidden_dim, 3)
    r1=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l2=Conv1D(hidden_dim, 3)
    r2=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l3=Conv1D(hidden_dim, 3)
    r3=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l4=Conv1D(hidden_dim, 3)
    r4=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l5=Conv1D(hidden_dim, 3)
    r5=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l6=Conv1D(hidden_dim, 3)
    s1=MaxPooling1D(3)(l1(seq_input1))
    s1=concatenate([r1(s1), s1])
    s1=MaxPooling1D(3)(l2(s1))
    s1=concatenate([r2(s1), s1])
    s1=MaxPooling1D(3)(l3(s1))
    s1=concatenate([r3(s1), s1])
    s1=MaxPooling1D(3)(l4(s1))
    s1=concatenate([r4(s1), s1])
    s1=MaxPooling1D(3)(l5(s1))
    s1=concatenate([r5(s1), s1])
    s1=l6(s1)
    s1=GlobalAveragePooling1D()(s1)
    s2=MaxPooling1D(3)(l1(seq_input2))
    s2=concatenate([r1(s2), s2])
    s2=MaxPooling1D(3)(l2(s2))
    s2=concatenate([r2(s2), s2])
    s2=MaxPooling1D(3)(l3(s2))
    s2=concatenate([r3(s2), s2])
    s2=MaxPooling1D(3)(l4(s2))
    s2=concatenate([r4(s2), s2])
    s2=MaxPooling1D(3)(l5(s2))
    s2=concatenate([r5(s2), s2])
    s2=l6(s2)
    s2=GlobalAveragePooling1D()(s2)
    merge_text = multiply([s1, s2])
    x = Dense(100, activation='linear')(merge_text)
    x = LeakyReLU(alpha=0.3)(x)
    x = Dense(int((hidden_dim+7)/2), activation='linear')(x)
    x = LeakyReLU(alpha=0.3)(x)
    main_output = Dense(2, activation='softmax')(x)
    merge_model = Model(inputs=[seq_input1, seq_input2], outputs=[main_output])
    return merge_model

def get_traintest_split(class_labels, train_length):

    train = class_labels[:train_length]
    test = class_labels[train_length:]
    train = []
    test = []
    for i in range(0, train_length):
        train.append(i)
    for j in range(train_length, len(class_labels)):
        test.append(j)

    return [(np.asarray(train), np.asarray(test))]

def get_crossvalidation_splits(class_labels, nsplits=5):
    kf = KFold(n_splits=nsplits, shuffle=True, random_state=10312020)
    tries = 5
    cur = 0
    train_test = []
    for train, test in kf.split(class_labels):
        if np.sum(class_labels[train], 0)[0] > 0.8 * len(train) or np.sum(class_labels[train], 0)[0] < 0.2 * len(train):
            continue
        train_test.append((train, test))
        cur += 1
        if cur >= tries:
            break
    return train_test


def get_test_results(id2index, raw_data, test_indices, class_labels, predictions):
    # id2index provides protein IDs from dictionary {id: number}
    # raw_data provides train+test data, interaction as protein numbers and label (0 or 1)
    # test_indices provides indices of test data within raw_data
    # class_labels provides labels for each interaction in raw_data encoded as [0, 1] or [1, 0]
    # predictions contains model predictions of each interaction for each label
    prob_results = []
    for ppi in range(0, len(test_indices)):
        proteinA_num = int(raw_data[test_indices[ppi]][0])
        proteinB_num = int(raw_data[test_indices[ppi]][1])
        #print('Label:', raw_data_test[ppi][-1])
        for k, v in id2index.items():
            if v == proteinA_num:
                proteinA = k
            if v == proteinB_num:
                proteinB = k  
        # Class label
        #true_interaction = class_labels[test_indices[ppi]][0]
        # Prediciton of positive interaction (index is 0 from CLASS_MAP)
        prob_pos_interaction = '{:f}'.format(float(predictions[ppi][0]))

        prob_results.append([proteinA + ' ' + proteinB + ' ' + str(prob_pos_interaction)])
        
    return np.asarray(prob_results, dtype=str)

def convert_num_to_protein(raw_data, id2index):
    # returns the raw_data with protein IDs instead of protein number from id2index
    converted = []
    for i in range(0, len(raw_data)):
        proteinA_num = raw_data[i][0]
        proteinB_num = raw_data[i][1]
        label = raw_data[i][2]
        for k, v in id2index.items():
            if v == proteinA_num:
                proteinA = k
            if v == proteinB_num:
                proteinB = k
        converted.append([proteinA, proteinB, label])
    return converted

if __name__ == "__main__":
    if args.cpu:
        print("\n---Troubleshooting with CPU only...model won't run---\n")
    else:
        # Setup GPU
        KTF.set_session(get_session())
    
    t_start = time()
    
    # Get protein sequences
    id2index = {}
    seqs = []
    index = 0
    for line in open(ID2SEQ_FILE):
        line = line.strip().split('\t')
        id2index[line[0]] = index
        seqs.append(line[1])
        index += 1
    print("Number of protein sequences:", index)
    
    sid1_index = 0
    sid2_index = 1
    
    if not CROSS_VALIDATE:
        # Process training data
        seq_array_train = []
        raw_data_train, seq_array_train, id2_aid_train = get_interaction_data(TRAIN_FILE, sid1_index, sid2_index, id2index, seqs, seq_array_train)
        seq_tensor_train = np.array([SEQ2T.embed_normalized(line, SEQ_SIZE) for line in tqdm(seq_array_train)])
        seq_index1_train = np.array([line[sid1_index] for line in tqdm(raw_data_train)])
        seq_index2_train = np.array([line[sid2_index] for line in tqdm(raw_data_train)])
        class_labels_train = create_class_labels(raw_data_train, CLASS_MAP, label_index)
        
        # Process testing data
        seq_array_test = []
        raw_data_test, seq_array_test, id2_aid_test = get_interaction_data(TEST_FILE, sid1_index, sid2_index, id2index, seqs, seq_array_test)
        seq_tensor_test = np.array([SEQ2T.embed_normalized(line, SEQ_SIZE) for line in tqdm(seq_array_test)])
        seq_index1_test = np.array([line[sid1_index] for line in tqdm(raw_data_test)])
        seq_index2_test = np.array([line[sid2_index] for line in tqdm(raw_data_test)])
        class_labels_test = create_class_labels(raw_data_test, CLASS_MAP, label_index)
        
        # Combine to common variable, but keep train/test split
        raw_data = np.concatenate((raw_data_train, raw_data_test))
        seq_tensor = np.concatenate((seq_tensor_train, seq_tensor_test))
        seq_index1 = np.concatenate((seq_index1_train, seq_index1_test))
        seq_index2 = np.concatenate((seq_index2_train, seq_index2_test))
        class_labels = np.concatenate((class_labels_train, class_labels_test))
        train_test = get_traintest_split(class_labels, class_labels_train.shape[0])
    else:
        seq_array = []
        raw_data, seq_array, id2_aid = get_interaction_data(TRAIN_FILE, sid1_index, sid2_index, id2index, seqs, seq_array)
        seq_tensor = np.array([SEQ2T.embed_normalized(line, SEQ_SIZE) for line in tqdm(seq_array)])
        seq_index1 = np.array([line[sid1_index] for line in tqdm(raw_data)])
        seq_index2 = np.array([line[sid2_index] for line in tqdm(raw_data)])
        class_labels = create_class_labels(raw_data, CLASS_MAP, label_index)
        train_test = get_crossvalidation_splits(class_labels, nsplits=5)

    if args.cpu:
        print("\nExiting before building model.")
        exit()
        
    avg_accuracy = []
    avg_precision = []
    avg_recall = []
    avg_specificity = []
    avg_f1 = []
    avg_mcc = []
    avg_roc_auc = []
    avg_pr_auc = []
        
    # Train and test model
    num_hit = num_total = num_pos = num_true_pos = num_false_pos = num_true_neg = num_false_neg = 0.
    batch_size1 = 256
    cv = 0
    for train, test in train_test:
        
        if not CROSS_VALIDATE and pretrained != None:
            merge_model = pickle.load(open(pretrained, 'rb'))
        else:
            merge_model = None
            merge_model = build_model()
            adam = Adam(lr=0.001, amsgrad=True, epsilon=1e-6)
            rms = RMSprop(lr=0.001)
            merge_model.compile(optimizer=rms, loss='categorical_crossentropy', metrics=['accuracy'])
            merge_model.fit([seq_tensor[seq_index1[train]], seq_tensor[seq_index2[train]]], class_labels[train], batch_size=batch_size1, epochs=n_epochs)
            
        if not CROSS_VALIDATE and args.saveModel:
            #pickle.dump(merge_model, open(os.getcwd()+'/Models/' + TRAIN_FILE.split('/')[-1].replace('.tsv', '_PIPR.model'), 'wb'))
            pickle.dump(merge_model, open(args.saveModel, "wb"))
            
        pred = merge_model.predict([seq_tensor[seq_index1[test]], seq_tensor[seq_index2[test]]])
        
        for i in range(len(class_labels[test])):        
            num_total += 1
            if np.argmax(class_labels[test][i]) == np.argmax(pred[i]):
                num_hit += 1
            if class_labels[test][i][0] > 0.:
                num_pos += 1.
                if pred[i][0] > pred[i][1]:
                    num_true_pos += 1
                else:
                    num_false_neg += 1
            else:
                if pred[i][0] > pred[i][1]:
                    num_false_pos += 1
                else:
                    num_true_neg += 1
        
        auc_roc_test = roc_auc_score(class_labels[test], pred)
        auc_pr_test = average_precision_score(class_labels[test], pred)
        
        if not CROSS_VALIDATE:
            # Save interaction probability results
            prob_results = get_test_results(id2_aid_test, raw_data, test, class_labels, pred)
            if args.predictions:
            #np.savetxt(os.getcwd()+'/Results/predictions_' + TRAIN_FILE.split('/')[-1].replace('.tsv', '_') + TEST_FILE.split('/')[-1].replace('.tsv', '.txt'), prob_results, fmt='%s', delimiter='\n'
                np.savetxt(args.predictions, prob_results, fmt='%s', delimiter='\n')
        else:
            # Save interaction probability results
            prob_results = get_test_results(id2_aid, raw_data, test, class_labels, pred)
            #np.savetxt(os.getcwd()+'/Results/predictions_' + TRAIN_FILE.split('/')[-1].replace('.tsv', '_') + TEST_FILE.split('/')[-1].replace('.tsv', '_') + 'fold-' + str(cv) + '.txt', prob_results, fmt='%s', delimiter='\n')
            if args.predictions:
                np.savetxt(args.predictions, prob_results, fmt='%s', delimiter='\n')
        
        print("======== Fold", cv)
        print('\ntp=%0.0f \nfp=%0.0f \ntn=%0.0f \nfn=%0.0f \n'%(num_true_pos, num_false_pos, num_true_neg, num_false_neg))
        cv += 1
        accuracy = (num_true_pos + num_true_neg) / num_total
        prec = num_true_pos / (num_true_pos + num_false_pos + 1e-06)
        recall = num_true_pos / (num_true_pos + num_false_neg + 1e-06)
        spec = num_true_neg / (num_true_neg + num_false_pos + 1e-06)
        f1 = 2. * (prec * recall) / (prec + recall + 1e-06)
        mcc = (num_true_pos * num_true_neg - num_false_pos * num_false_neg) / (((num_true_pos + num_false_pos + 1e-06) * (num_true_pos + num_false_neg + 1e-06) * (num_false_pos + num_true_neg + 1e-06) * (num_true_neg + num_false_neg + 1e-06)) ** 0.5)
        print('acc=', accuracy, '\nprec=', prec, '\nrecall=', recall, '\nspec=', spec, '\nf1=', f1, '\nmcc=', mcc)
        print('auc_roc=', auc_roc_test, '\nauc_pr=', auc_pr_test)
        
        avg_accuracy.append(accuracy)
        avg_precision.append(prec)
        avg_recall.append(recall)
        avg_specificity.append(spec)
        avg_f1.append(f1)
        avg_mcc.append(mcc)
        avg_roc_auc.append(auc_roc_test)
        avg_pr_auc.append(auc_pr_test)
        print('\n', time() - t_start, 'seconds to complete')
    
    # Write results to file
    if rst_file:
        with open(rst_file, 'w') as fp:
            fp.write(('accuracy=%.4f (+/- %.4f)'%(np.mean(avg_accuracy), np.std(avg_accuracy))
                      + '\nprecision=%.4f (+/- %.4f)'%(np.mean(avg_precision), np.std(avg_precision)) 
                      + '\nrecall=%.4f (+/- %.4f)'%(np.mean(avg_recall), np.std(avg_recall)) 
                      + '\nspecificity=%.4f (+/- %.4f)'%(np.mean(avg_specificity), np.std(avg_specificity)) 
                      + '\nf1=%.4f (+/- %.4f)'%(np.mean(avg_f1), np.std(avg_f1)) 
                      + '\nmcc=%.4f (+/- %.4f)'%(np.mean(avg_mcc), np.std(avg_mcc))
                      + '\nroc_auc=%.4f (+/- %.4f)' % (np.mean(avg_roc_auc), np.std(avg_roc_auc))
                      + '\npr_auc=%.4f (+/- %.4f)' % (np.mean(avg_pr_auc), np.std(avg_pr_auc))
                      + '\ntime=%.2f'%(time()-t_start)
                      + '\n'))
