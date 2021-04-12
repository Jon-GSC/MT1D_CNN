#load functions for cnn training and plot.
#
import os
import numpy as np
import six

import matplotlib.pyplot as plt
import matplotlib.tri as tri
plt.style.use('seaborn-white')
import seaborn as sns
sns.set_style("white")

from scipy.interpolate import interp1d, interp2d, splev, splrep, Rbf, InterpolatedUnivariateSpline
from scipy.spatial import Delaunay
import keras
from keras.layers import Input,Dropout,BatchNormalization,Activation,Add,Lambda,Flatten,Dense,Reshape
from keras import backend as K
from keras.regularizers import l2
from keras import optimizers
from keras.engine.topology import Input
from keras.engine.training import Model
from keras.layers.convolutional import Conv1D, UpSampling1D, Conv2DTranspose
from keras.layers.core import Activation, SpatialDropout2D
from keras.layers.merge import concatenate,add
from keras.layers.normalization import BatchNormalization
from keras.layers.pooling import MaxPooling1D
from keras.losses import binary_crossentropy

import tensorflow as tf

CHANNEL_AXIS = 2
input_size = 256
#-----------------------------------------------------------------------------------------------------------------------
#Params and helpers
def exp10(my_list):  #
    return [ 10**x for x in my_list ]

def Data_Augumentation(Depth, x_train, y_train):  #
    d_size = x_train.shape

    x_train_ex = np.empty([0,input_size])
    y_train_ex = np.empty([0,input_size])
    for i in range(d_size[0]):
        Dxmaxmin = Depth[x_train[i].reshape(-1) > 0.0]
        Dymaxmin = Depth[y_train[i].reshape(-1) > 0.0]
        k = 0
        D_start = min(Dxmaxmin)
        temp0 = x_train[i].reshape(-1)
        temp2 = y_train[i].reshape(-1)
        while D_start < min(Dymaxmin) and k < 6:
            temp1 = np.concatenate((temp0[100:],temp0[0:100]))
            temp0 = temp1.copy()
            x_train_ex = np.vstack((x_train_ex, temp1))  # add to x_train.
            temp3 = np.concatenate((temp2[100:],temp2[0:100]))
            temp2 = temp3.copy()
            y_train_ex = np.vstack((y_train_ex,temp3))  # add to y_train.

            D_start += 20
            k += 1
            pass
    x_train = np.concatenate((x_train.reshape(-1,input_size), x_train_ex), axis=0).reshape(-1,input_size,1)
    y_train = np.concatenate((y_train.reshape(-1,input_size), y_train_ex), axis=0).reshape(-1,input_size,1)
    return x_train, y_train

def BatchActivate(x):
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    return x

def convolution_block(x, filters, size, strides=1, padding='same', activation=True):
    x = Conv1D(filters, size, strides=strides, padding=padding)(x)
    if activation == True:
        x = BatchActivate(x)
    return x


def residual_block(blockInput, num_filters=16, batch_activate=False):
    x = BatchActivate(blockInput)
    x = convolution_block(x, num_filters, 3)
    x = convolution_block(x, num_filters, 3, activation=False)
    x = Add()([x, blockInput])
    if batch_activate:
        x = BatchActivate(x)
    return x


def Conv1DTranspose(input_tensor, filters, ksize, strides=2,padding='same'):
    x = Lambda(lambda x: K.expand_dims(x,axis=2))(input_tensor)
    x = Conv2DTranspose(filters=filters,kernel_size=(ksize,1),strides=(strides,1),padding=padding)(x)
    x = Lambda(lambda x: K.squeeze(x, axis=2))(x)
    return x

def shape_match(input_tensor,start_neurons,size):
    return Conv1D(start_neurons, size, activation=None, padding="valid")(input_tensor)

def build_model(input_layer, output_size, start_neurons, size=3, DropoutRatio=0.5):
    conv1 = Conv1D(start_neurons * 1, size, activation="relu", padding="same")(input_layer)
    conv1 = residual_block(conv1, start_neurons * 1)
    conv1 = residual_block(conv1, start_neurons * 1, True)
    pool1 = MaxPooling1D(2)(conv1)
    pool1 = Dropout(DropoutRatio / 2)(pool1)

    conv2 = Conv1D(start_neurons * 2, size, activation=None, padding="same")(pool1)
    conv2 = residual_block(conv2, start_neurons * 2)
    conv2 = residual_block(conv2, start_neurons * 2, True)
    pool2 = MaxPooling1D(2)(conv2)
    pool2 = Dropout(DropoutRatio)(pool2)

    conv3 = Conv1D(start_neurons * 4, size, activation=None, padding="same")(pool2)
    conv3 = residual_block(conv3, start_neurons * 4)
    conv3 = residual_block(conv3, start_neurons * 4, True)
    pool3 = MaxPooling1D(2)(conv3)
    pool3 = Dropout(DropoutRatio)(pool3)

    conv4 = Conv1D(start_neurons * 8, size, activation=None, padding="same")(pool3)
    conv4 = residual_block(conv4, start_neurons * 8)
    conv4 = residual_block(conv4, start_neurons * 8, True)
    pool4 = MaxPooling1D(2)(conv4)
    pool4 = Dropout(DropoutRatio)(pool4)

    convm = Conv1D(start_neurons * 16, size, activation=None, padding="same")(pool4)
    convm = residual_block(convm, start_neurons * 16)
    convm = residual_block(convm, start_neurons * 16, True)

    deconv4 = Conv1DTranspose(convm, start_neurons * 8, size, strides=2, padding="same")
    uconv4 = concatenate([deconv4, conv4])
    uconv4 = Dropout(DropoutRatio)(uconv4)

    uconv4 = Conv1D(start_neurons * 8, size, activation=None, padding="same")(uconv4)
    uconv4 = residual_block(uconv4, start_neurons * 8)
    uconv4 = residual_block(uconv4, start_neurons * 8, True)

    deconv3 = Conv1DTranspose(uconv4,start_neurons * 4, size, strides=2, padding="same")
    uconv3 = concatenate([deconv3, conv3])
    uconv3 = Dropout(DropoutRatio)(uconv3)

    uconv3 = Conv1D(start_neurons * 4, size, activation=None, padding="same")(uconv3)
    uconv3 = residual_block(uconv3, start_neurons * 4)
    uconv3 = residual_block(uconv3, start_neurons * 4, True)

    deconv2 = Conv1DTranspose(uconv3, start_neurons * 2, size, strides=2, padding="same")
    uconv2 = concatenate([deconv2, conv2])

    uconv2 = Dropout(DropoutRatio)(uconv2)
    uconv2 = Conv1D(start_neurons * 2, size, activation=None, padding="same")(uconv2)
    uconv2 = residual_block(uconv2, start_neurons * 2)
    uconv2 = residual_block(uconv2, start_neurons * 2, True)

    deconv1 = Conv1DTranspose(uconv2, start_neurons * 1, size, strides=2, padding="same")
    uconv1 = concatenate([deconv1, conv1])

    uconv1 = Dropout(DropoutRatio)(uconv1)
    uconv1 = Conv1D(start_neurons * 1, size, activation=None, padding="same")(uconv1)
    uconv1 = residual_block(uconv1, start_neurons * 1)
    uconv1 = residual_block(uconv1, start_neurons * 1, True)

    output_layer_noActi = Conv1D(1, 1, padding="same", activation=None)(uconv1)
    uconv0 = Flatten()(output_layer_noActi)
    uconv0 = Dense(output_size, activation=None)(uconv0)
    model = Model(input_layer, uconv0)
    return model

# Get Iou Vector
def get_iou_vector(A, B):
    batch_size = A.shape[0]
    metric = []
    for batch in range(batch_size):
        t, p = A[batch] > 0, B[batch] > 0
        intersection = np.logical_and(t, p)
        union = np.logical_or(t, p)
        iou = (np.sum(intersection > 0) + 1e-15) / (np.sum(union > 0) + 1e-15)
        thresholds = np.arange(0.0001, 0.2, 0.0001)    # how many values and range is important....
        s = []
        for thresh in thresholds:
            s.append(iou > thresh)
        metric.append(np.mean(s))
    return np.mean(metric)

def my_iou_metric(label, pred):
    return tf.compat.v1.py_func(get_iou_vector, [label, pred > 0.05], tf.float64)

# loss functions
def dice_coef(y_true, y_pred):
    y_true_f = K.flatten(y_true)
    y_pred = K.cast(y_pred, 'float32')
    y_pred_f = K.cast(K.greater(K.flatten(y_pred), 0.5), 'float32')
    intersection = y_true_f * y_pred_f
    score = 2. * K.sum(intersection) / (K.sum(y_true_f) + K.sum(y_pred_f))
    return score

# https://gist.github.com/bshishov/5dc237f59f019b26145648e2124ca1c9 .revised by Jon,
def _error(y_true, y_pred):
    """ Simple error """
    return y_true - y_pred

def mse(y_true, y_pred):
    """ Mean Squared Error """
    return K.mean(K.square(_error(y_true, y_pred)))

def rmse(y_true, y_pred):
    """ Root Mean Squared Error """
    return K.sqrt(mse(y_true, y_pred))

def nrmse(y_true, y_pred):
    """ Normalized Root Mean Squared Error """
    # normalized by corelation?
    return rmse(y_true, y_pred) / (K.max(y_true) - K.min(y_true)) #/ corrcoef0   #??

def rmse_m(y_true, y_pred):
    """ Normalized Root Mean Squared Error plus minimum model """
    return rmse(y_true, y_pred) + rmse(y_pred,0.0)/50/10   # add model-minimum,

def root_mean_squared_error(y_true, y_pred):
    return K.sqrt(K.mean(K.square(y_pred - y_true)))

def correlation_coefficient_loss(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x)
    my = K.mean(y)
    xm, ym = x-mx, y-my
    r_num = K.sum(tf.multiply(xm,ym))
    r_den = K.sqrt(tf.multiply(K.sum(K.square(xm)), K.sum(K.square(ym))))
    r = r_num / r_den
    r = K.maximum(K.minimum(r, 1.0), -1.0)
    return 1 - K.square(r)

def dice_loss(y_true, y_pred):
    smooth = 1.
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = y_true_f * y_pred_f
    score = (2. * K.sum(intersection) + smooth) / (K.sum(y_true_f) + K.sum(y_pred_f) + smooth)
    return 1. - score


def bce_dice_loss(y_true, y_pred):
    return binary_crossentropy(y_true, y_pred) + dice_loss(y_true, y_pred)


def bce_logdice_loss(y_true, y_pred):
    return binary_crossentropy(y_true, y_pred) - K.log(1. - dice_loss(y_true, y_pred))


def weighted_bce_loss(y_true, y_pred, weight):
    epsilon = 1e-7
    y_pred = K.clip(y_pred, epsilon, 1. - epsilon)
    logit_y_pred = K.log(y_pred / (1. - y_pred))
    loss = weight * (logit_y_pred * (1. - y_true) +
                     K.log(1. + K.exp(-K.abs(logit_y_pred))) + K.maximum(-logit_y_pred, 0.))
    return K.sum(loss) / K.sum(weight)


def weighted_dice_loss(y_true, y_pred, weight):
    smooth = 1.
    w, m1, m2 = weight, y_true, y_pred
    intersection = (m1 * m2)
    score = (2. * K.sum(w * intersection) + smooth) / (K.sum(w * m1) + K.sum(w * m2) + smooth)
    loss = 1. - K.sum(score)
    return loss


# code download from: https://github.com/bermanmaxim/LovaszSoftmax
def lovasz_grad(gt_sorted):
    """
    Computes gradient of the Lovasz extension w.r.t sorted errors
    See Alg. 1 in paper
    """
    gts = tf.reduce_sum(gt_sorted)
    intersection = gts - tf.cumsum(gt_sorted)
    union = gts + tf.cumsum(1. - gt_sorted)
    jaccard = 1. - intersection / union
    jaccard = tf.concat((jaccard[0:1], jaccard[1:] - jaccard[:-1]), 0)
    return jaccard


# --------------------------- BINARY LOSSES ---------------------------
def lovasz_hinge(logits, labels, per_image=True, ignore=None):
    """
    Binary Lovasz hinge loss
      logits: [B, H, W] Variable, logits at each pixel (between -\infty and +\infty)
      labels: [B, H, W] Tensor, binary ground truth masks (0 or 1)
      per_image: compute the loss per image instead of per batch
      ignore: void class id
    """
    if per_image:
        def treat_image(log_lab):
            log, lab = log_lab
            log, lab = tf.expand_dims(log, 0), tf.expand_dims(lab, 0)
            log, lab = flatten_binary_scores(log, lab, ignore)
            return lovasz_hinge_flat(log, lab)

        losses = tf.map_fn(treat_image, (logits, labels), dtype=tf.float32)
        loss = tf.reduce_mean(losses)
    else:
        loss = lovasz_hinge_flat(*flatten_binary_scores(logits, labels, ignore))
    return loss


def lovasz_hinge_flat(logits, labels):
    """
    Binary Lovasz hinge loss
      logits: [P] Variable, logits at each prediction (between -\infty and +\infty)
      labels: [P] Tensor, binary ground truth labels (0 or 1)
      ignore: label to ignore
    """
    def compute_loss():
        labelsf = tf.cast(labels, logits.dtype)
        signs = 2. * labelsf - 1.
        errors = 1. - logits * tf.stop_gradient(signs)
        errors_sorted, perm = tf.nn.top_k(errors, k=tf.shape(errors)[0], name="descending_sort")
        gt_sorted = tf.gather(labelsf, perm)
        grad = lovasz_grad(gt_sorted)
        loss = tf.tensordot(tf.nn.elu(errors_sorted) + 1, tf.stop_gradient(grad), 1, name="loss_non_void")
        return loss

    # deal with the void prediction case (only void pixels)
    loss = tf.cond(tf.equal(tf.shape(logits)[0], 0),
                   lambda: tf.reduce_sum(logits) * 0.,
                   compute_loss,
                   strict=True,
                   name="loss"
                   )
    return loss

def flatten_binary_scores(scores, labels, ignore=None):
    """
    Flattens predictions in the batch (binary case)
    Remove labels equal to 'ignore'
    """
    scores = tf.reshape(scores, (-1,))
    labels = tf.reshape(labels, (-1,))
    if ignore is None:
        return scores, labels
    valid = tf.not_equal(labels, ignore)
    vscores = tf.boolean_mask(scores, valid, name='valid_scores')
    vlabels = tf.boolean_mask(labels, valid, name='valid_labels')
    return vscores, vlabels


def lovasz_loss(y_true, y_pred):
    y_true, y_pred = K.cast(K.squeeze(y_true, -1), 'int32'), K.cast(K.squeeze(y_pred, -1), 'float32')
    # logits = K.log(y_pred / (1. - y_pred))
    logits = y_pred  # Jiaxin
    loss = lovasz_hinge(logits, y_true, per_image=True, ignore=None)
    return loss


# ResNet 34
# https://github.com/raghakot/keras-resnet/blob/master/resnet.py
def _bn_relu(input):
    """Helper to build a BN -> relu block
    """
    norm = BatchNormalization(axis=2)(input)
    return Activation("relu")(norm)


def _conv_bn_relu(**conv_params):
    """Helper to build a conv -> BN -> relu block
    """
    filters = conv_params["filters"]
    kernel_size = conv_params["kernel_size"]
    strides = conv_params.setdefault("strides", 1)
    kernel_initializer = conv_params.setdefault("kernel_initializer", "he_normal")
    padding = conv_params.setdefault("padding", "same")
    kernel_regularizer = conv_params.setdefault("kernel_regularizer", l2(1.e-4))

    def f(input):
        conv = Conv1D(filters=filters, kernel_size=kernel_size,
                      strides=strides, padding=padding,
                      kernel_initializer=kernel_initializer,
                      kernel_regularizer=kernel_regularizer)(input)
        return _bn_relu(conv)
    return f


def _bn_relu_conv(**conv_params):
    """Helper to build a BN -> relu -> conv block.
    This is an improved scheme proposed in http://arxiv.org/pdf/1603.05027v2.pdf
    """
    filters = conv_params["filters"]
    kernel_size = conv_params["kernel_size"]
    strides = conv_params.setdefault("strides", 1)
    kernel_initializer = conv_params.setdefault("kernel_initializer", "he_normal")
    padding = conv_params.setdefault("padding", "same")
    kernel_regularizer = conv_params.setdefault("kernel_regularizer", l2(1.e-4))

    def f(input):
        activation = _bn_relu(input)
        return Conv1D(filters=filters, kernel_size=kernel_size,
                      strides=strides, padding=padding,
                      kernel_initializer=kernel_initializer,
                      kernel_regularizer=kernel_regularizer)(activation)
    return f


def _shortcut(input, residual):
    """Adds a shortcut between input and residual block and merges them with "sum"
    """
    # Expand channels of shortcut to match residual.
    # Stride appropriately to match residual (width, height)
    # Should be int if network architecture is correctly configured.
    CHANNEL_AXIS = 2

    input_shape = K.int_shape(input)
    residual_shape = K.int_shape(residual)
    stride_width = int(round(input_shape[ROW_AXIS] / residual_shape[ROW_AXIS]))
    stride_height = int(round(input_shape[COL_AXIS] / residual_shape[COL_AXIS]))
    equal_channels = input_shape[CHANNEL_AXIS] == residual_shape[CHANNEL_AXIS]

    shortcut = input
    # 1 X 1 conv if shape is different. Else identity.
    if stride_width > 1 or not equal_channels:
        shortcut = Conv1D(filters=residual_shape[CHANNEL_AXIS],
                          kernel_size=1,
                          strides=stride_width,
                          padding="valid",
                          kernel_initializer="he_normal",
                          kernel_regularizer=l2(0.0001))(input)
    return add([shortcut, residual])


def basic_block(filters, init_strides=1, is_first_block_of_first_layer=False):
    """Basic 3 X 3 convolution blocks for use on resnets with layers <= 34.
    """
    def f(input):

        if is_first_block_of_first_layer:
            # don't repeat bn->relu since we just did bn->relu->maxpool
            conv1 = Conv1D(filters=filters, kernel_size=3,
                           strides=init_strides,
                           padding="same",
                           kernel_initializer="he_normal",
                           kernel_regularizer=l2(1e-4))(input)
        else:
            conv1 = _bn_relu_conv(filters=filters, kernel_size=3,
                                  strides=init_strides)(input)

        residual = _bn_relu_conv(filters=filters, kernel_size=3)(conv1)
        return _shortcut(input, residual)
    return f


def _residual_block(block_function, filters, repetitions, is_first_layer=False):
    """Builds a residual block with repeating bottleneck blocks.
    """
    def f(input):
        for i in range(repetitions):
            init_strides = 1
            if i == 0 and not is_first_layer:
                init_strides = 2
            input = block_function(filters=filters, init_strides=init_strides,
                                   is_first_block_of_first_layer=(is_first_layer and i == 0))(input)
        return input
    return f


def _handle_dim_ordering():
    global ROW_AXIS
    global COL_AXIS
    global CHANNEL_AXIS
    if K.common.image_dim_ordering() == 'tf':
        ROW_AXIS = 1
        COL_AXIS = 2
        CHANNEL_AXIS = 3
    else:
        CHANNEL_AXIS = 1
        ROW_AXIS = 2
        COL_AXIS = 3


def _get_block(identifier):
    if isinstance(identifier, six.string_types):
        res = globals().get(identifier)
        if not res:
            raise ValueError('Invalid {}'.format(identifier))
        return res
    return identifier


class ResnetBuilder(object):
    @staticmethod
    def build(input_shape, block_fn, repetitions, input_tensor):
        _handle_dim_ordering()
        if len(input_shape) != 2:
            raise Exception("Input shape should be a tuple (nb_channels, nb_rows, nb_cols)")

        # Permute dimension order if necessary
        if K.common.image_dim_ordering() == 'tf':
            input_shape = (input_shape[1], input_shape[0])

        # Load function from str if needed.
        block_fn = _get_block(block_fn)

        if input_tensor is None:
            img_input = Input(shape=input_shape)
        else:
            if not K.is_keras_tensor(input_tensor):
                img_input = Input(tensor=input_tensor, shape=input_shape)
            else:
                img_input = input_tensor

        conv1 = _conv_bn_relu(filters=64, kernel_size=7, strides=2)(img_input)
        pool1 = MaxPooling1D(pool_size=3, strides=2, padding="same")(conv1)

        block = pool1
        filters = 64
        for i, r in enumerate(repetitions):
            block = _residual_block(block_fn, filters=filters, repetitions=r, is_first_layer=(i == 0))(block)
            filters *= 2
        # Last activation
        block = _bn_relu(block)
        model = Model(inputs=img_input, outputs=block)
        return model
    @staticmethod
    def build_resnet_34(input_shape, input_tensor):
        return ResnetBuilder.build(input_shape, basic_block, [3, 4, 6, 3], input_tensor)


def UResNet34(input_shape=None, output_size=None, classes=1, decoder_filters=32, input_tensor=None, activation='sigmoid', **kwargs):
    input_layer0 = Input(input_shape)
    input_layer1 = Input(input_shape)

    output_layer0 = build_model(input_layer0, output_size, 32, 0.5)
    output_layer1 = build_model(input_layer1, output_size, 32, 0.5)

    model0 = Model(input_layer0, output_layer0)
    model1 = Model(input_layer1, output_layer1)

    combined = concatenate([model0.output, model1.output])
    uconv0 = Dense(output_size, activation=None, name='predictions')(combined)
    uconv0 = Reshape((uconv0.shape[1],1))(uconv0)
    model = Model([model0.input,model1.input], uconv0)
    model.name = 'unet-combined-Jon'
    return model

"""
used for converting the decoded image to rle mask Fast compared to previous one
"""
def rle_encode(im):
    '''
    im: numpy array, 1 - mask, 0 - background
    Returns run length as string formated
    '''
    pixels = im.flatten(order = 'F')
    pixels = np.concatenate([[0], pixels, [0]])
    runs = np.where(pixels[1:] != pixels[:-1])[0] + 1
    runs[1::2] -= runs[::2]
    return ' '.join(str(x) for x in runs)


# Score the model and do a threshold optimization by the best IoU.
# src: https://www.kaggle.com/aglotero/another-iou-metric
def iou_metric(y_true_in, y_pred_in, print_table=False):
    labels = y_true_in
    y_pred = y_pred_in

    true_objects = 2
    pred_objects = 2

    temp1 = np.histogram2d(labels.flatten(), y_pred.flatten(), bins=([0, 0.5, 1], [0, 0.5, 1]))
    intersection = temp1[0]
    area_true = np.histogram(labels, bins=[0, 0.5, 1])[0]
    area_pred = np.histogram(y_pred, bins=[0, 0.5, 1])[0]
    area_true = np.expand_dims(area_true, -1)
    area_pred = np.expand_dims(area_pred, 0)
    # Compute union
    union = area_true + area_pred - intersection
    # Exclude background from the analysis
    intersection = intersection[1:, 1:]
    intersection[intersection == 0] = 1e-9
    union = union[1:, 1:]
    union[union == 0] = 1e-9

    # Compute the intersection over union
    iou = intersection / union

    # Precision helper function
    def precision_at(threshold, iou):
        matches = iou > threshold
        true_positives = np.sum(matches, axis=1) == 1  # Correct objects
        false_positives = np.sum(matches, axis=0) == 0  # Missed objects
        false_negatives = np.sum(matches, axis=1) == 0  # Extra objects
        tp, fp, fn = np.sum(true_positives), np.sum(false_positives), np.sum(false_negatives)
        return tp, fp, fn

    # Loop over IoU thresholds
    prec = []
    if print_table:
        print("Thresh\tTP\tFP\tFN\tPrec.")
    for t in np.arange(0.5, 1.0, 0.05):
        tp, fp, fn = precision_at(t, iou)
        if (tp + fp + fn) > 0:
            p = tp / (tp + fp + fn)
        else:
            p = 0
        if print_table:
            print("{:1.3f}\t{}\t{}\t{}\t{:1.3f}".format(t, tp, fp, fn, p))
        prec.append(p)

    if print_table:
        print("AP\t-\t-\t-\t{:1.3f}".format(np.mean(prec)))
    return np.mean(prec)


def iou_metric_batch(y_true_in, y_pred_in):
    batch_size = y_true_in.shape[0]
    metric = []
    for batch in range(batch_size):
        value = iou_metric(y_true_in[batch], y_pred_in[batch])
        metric.append(value)
    return np.mean(metric)


def predict_result(model,x_test,img_size_target):
    preds_test = model.predict(x_test).reshape(-1, img_size_target)
    return preds_test

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'
#-----------------------------------------------------------------------------------------------------------------------
#
def remove_files(images_path_tmp):
    for root, dirs, files in os.walk(images_path_tmp, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))

def read_startup_0(startup_fn):   # edited from mtpy
    """
    reads in a 1D input file
    Arguments:
    ---------
        **inputfn** : full path to input file
    Returns:
    --------
        **Occam1D.indict** : dictionary with keys following the header and
            *'res'* : an array of resistivity values
    """
    if startup_fn is not None:
        startup_fn = startup_fn
    if startup_fn is None:
        raise IOError('Need to input a startup file.')

    startup_bn = os.path.basename(startup_fn)
    save_path = os.path.dirname(startup_fn)

    infid = open(startup_fn, 'r')
    ilines = infid.readlines()
    infid.close()
    indict = {}
    res = []
    # split the keys and values from the header information
    for iline in ilines:
        if iline.find(':') >= 0:
            ikey = iline[0:20].strip()[:-1]
            ivalue = iline[20:].split('!')[0].strip()
            if ikey.find('!') == 0:
                pass
            else:
                ikey.lower().replace(' ', '_')
                # setattr(self, ikey.lower().replace(' ', '_'), ivalue)
            indict[ikey[:-1]] = ivalue
        else:
            try:
                #print([float(xx) for xx in iline.split()])
                [res.append(float(xx)/1.0) for xx in iline.split()]
            except ValueError:
                pass
    return res

def transfer2occam_jon(resp_fn):  # revised from mtpy
    mode = 'det'
    # read original..
    if resp_fn is None:
        raise IOError('Need to input response file')
    # --> read response file
    dfid = open(resp_fn, 'r')
    dlines0 = dfid.readlines()
    dfid.close()
    finddict = {}
    for ii, dline in enumerate(dlines0):
        if dline.find('#') <= 3:
            fkey = dline[2:].strip().split(':')[0]
            fvalue = ii
            finddict[fkey] = fvalue

    # --> write file as a list of lines
    dlines = []
    dlines.append('Format:  EMData_1.1 \n')
    dlines.append('!mode:   {0}\n'.format(mode.upper()))
    dlines.append('!rotation_angle = {0:.2f}\n'.format(0.0))
    # needs a transmitter to work so put in a dummy one
    dlines.append('# Transmitters: 1\n')
    dlines.append('0 0 0 0 0 \n')

    freqency_lines = dlines0[finddict['Frequencies'] + 1:finddict['Receivers']]
    nf = len(freqency_lines)
    # write frequencies
    dlines.append('# Frequencies:   {0}\n'.format(nf))

    for fre_temp in freqency_lines:
        dlines.append(fre_temp)   # append frequency

    # needs a receiver to work so put in a dummy one
    dlines.append('# Receivers: 1 \n')
    dlines.append('0 0 0 0 0 0 \n')
    # write data
    dlines.append('# Data:{0}\n'.format(2 * nf))
    dlines.append('!{0}\n'.format('     '.join(['Type', 'Freq#', 'TX#', 'Rx#', 'Data', 'Std_Error'])))
    for dline in dlines0[finddict['Data'] + 1:]:
        if dline.find('!') == 0:
            pass
        else:
            dlst = dline.strip().split()
            if len(dlst) > 4:
                del dlst[4:6]
            dlines.append('    '.join(dlst)+'\n')
    # --> write file
    with open(resp_fn, 'w') as dfid:
        dfid.writelines(dlines)
    print('Transfer Data File for inversion : {0} by Jon'.format(resp_fn))


def transfer2occam_jon_1(input_size,resp_path0,savepath1,appres,phase):
    mode = 'det'
    fn = '/Occam1d_DataFile_DET.dat'
    resp_fn = resp_path0 + fn
    resp_fn_save = savepath1 + fn
# read original..
    if resp_fn is None:
        raise IOError('Need to input response file')
    # --> read response file
    dfid = open(resp_fn, 'r')
    dlines0 = dfid.readlines()
    dfid.close()
    finddict = {}
    for ii, dline in enumerate(dlines0):
        if dline.find('#') <= 3:
            fkey = dline[2:].strip().split(':')[0]
            fvalue = ii
            finddict[fkey] = fvalue
    # --> write file as a list of lines
    dlines = []
    dlines.append('Format:  EMData_1.1 \n')
    dlines.append('!mode:   {0}\n'.format(mode.upper()))
    dlines.append('!rotation_angle = {0:.2f}\n'.format(0.0))
    # needs a transmitter to work so put in a dummy one
    dlines.append('# Transmitters: 1\n')
    dlines.append('0 0 0 0 0 \n')
    freqency_lines = dlines0[finddict['Frequencies'] + 1:finddict['Receivers']]
    freq0 = np.logspace(np.log10(float(freqency_lines[0])), np.log10(float(freqency_lines[-1])), input_size).tolist()
    nf = len(freq0)
    dlines.append('# Frequencies:   {0}\n'.format(nf))
    for fre_temp in freq0:
        dlines.append(str(fre_temp)+'\n')   # append frequency
    # needs a receiver to work so put in a dummy one
    dlines.append('# Receivers: 1 \n')
    dlines.append('0 0 0 0 0 0 \n')
    # write data
    dlines.append('# Data:{0}\n'.format(2 * nf))
    dlines.append('!{0}\n'.format('     '.join(['Type', 'Freq#', 'TX#', 'Rx#', 'Data', 'Std_Error'])))
    for dline in dlines0[finddict['Data'] + 1:]:
        if dline.find('!') == 0:
            pass
        else:
            for idx0 in range(nf):
                dlst = ['103', str(idx0+1), str(0), str(1), str(appres[idx0]), str(1)]
                dlines.append('    '.join(dlst) + '\n')
                dlst = ['104', str(idx0+1), str(0), str(1), str(phase[idx0]), str(1)]
                dlines.append('    '.join(dlst) + '\n')
            break
    # --> write file
    with open(resp_fn_save, 'w') as dfid:
        dfid.writelines(dlines)
    print('Transfer Data File for inversion : {0} by Jon'.format(resp_fn_save))


def read_resp_file_0(input_size=128,resp_fn=None):      # edited from mtpy
    """
    read response file
     Arguments:
    ---------
        **resp_fn** : full path to response file
         **data_fn** : full path to data file
     Fills:
    --------

        *freq* : an array of frequencies with length nf

        *res_te* : TE resistivity array with shape (nf,4) for (0) data,
                  (1) dataerr, (2) model, (3) modelerr
         *res_tm* : TM resistivity array with shape (nf,4) for (0) data,
                  (1) dataerr, (2) model, (3) modelerr
        *phase_te* : TE phase array with shape (nf,4) for (0) data,
                    (1) dataerr, (2) model, (3) modelerr
        *phase_tm* : TM phase array with shape (nf,4) for (0) data,
                    (1) dataerr, (2) model, (3) modelerr
    """
    mode = 'det'
    res_te = []
    phase_te = []

    if resp_fn is None:
        raise IOError('Need to input response file')
    # --> read response file
    dfid = open(resp_fn, 'r')

    dlines = dfid.readlines()
    dfid.close()

    finddict = {}
    for ii, dline in enumerate(dlines):
        if dline.find('#') <= 3:
            fkey = dline[2:].strip().split(':')[0]
            fvalue = ii
            finddict[fkey] = fvalue

    freq = np.array([float(ff) for ff in dlines[finddict['Frequencies'] + 1:finddict['Receivers']]])

    for dline in dlines[finddict['Data'] + 1:]:
        if dline.find('!') == 0:
            pass
        else:
            dlst = dline.strip().split()
            if len(dlst) > 4:
                jj = int(dlst[1]) - 1
                dvalue = float(dlst[4])
                derr = float(dlst[5])
                rvalue = float(dlst[6])
                try:
                    rerr = float(dlst[7])
                except ValueError:
                    rerr = 1000.
                if dlst[0] == 'RhoZxy' or dlst[0] == '103':
                    res_te.append(rvalue)
                if dlst[0] == 'PhsZxy' or dlst[0] == '104':
                    phase_te.append(rvalue)

    freq0 = np.logspace(np.log10(freq[0]),np.log10(freq[-1]), input_size)   # interpolation with frequency,.

    f1d1 = interp1d(freq, np.array(res_te), kind='linear',fill_value='extrapolate')
    res_te = f1d1(freq0)
    f1d2 = interp1d(freq, np.array(phase_te), kind='linear',fill_value='extrapolate')
    phase_te = f1d2(freq0)

    appres_i = res_te.tolist()  # .
    phase_i = phase_te.tolist()  # .
    return appres_i,phase_i   #return app-res and phase

def calculate_misfit(input_size=128,resp_fn=None):
    mode = 'det'
    res_te = []
    phase_te = []

    if resp_fn is None:
        raise IOError('Need to input response file')
    # --> read response file
    dfid = open(resp_fn, 'r')

    dlines = dfid.readlines()
    dfid.close()

    finddict = {}
    for ii, dline in enumerate(dlines):
        if dline.find('#') <= 3:
            fkey = dline[2:].strip().split(':')[0]
            fvalue = ii
            finddict[fkey] = fvalue

    freq = np.array([float(ff) for ff in dlines[finddict['Frequencies'] + 1:finddict['Receivers']]])   #can optimize

    for dline in dlines[finddict['Data'] + 1:]:
        if dline.find('!') == 0:
            pass
        else:
            dlst = dline.strip().split()
            if len(dlst) > 4:
                jj = int(dlst[1]) - 1
                dvalue = float(dlst[4])
                derr = float(dlst[5])
                rvalue = float(dlst[6])
                try:
                    rerr = float(dlst[7])
                except ValueError:
                    rerr = 1000.
                if dlst[0] == 'RhoZxy' or dlst[0] == '103' or dlst[0] == '105':
                    res_te.append(abs(np.log10(dvalue)-np.log10(rvalue)))
                if dlst[0] == 'PhsZxy' or dlst[0] == '104' or dlst[0] == '106':
                    phase_te.append(abs(dvalue-rvalue)/90)

    misfit = np.mean(np.square(res_te))+np.mean(np.square(phase_te))
    return misfit


def res_generator_0(n_layer,ns=10):   #generate resistivity model
    nseed = np.random.randint(888)
    np.random.seed(nseed)
    icase = np.random.randint(2,4)
    if icase==1:
        res_i = 2.5 + 2.25*np.random.randn(n_layer-1)
        res_i[res_i < 0.0] = 0.0;    res_i[res_i > 6.0] = 6.0
        res_i = res_i - min(res_i)
    elif icase==2 or icase==3:
        t0 = np.linspace(0,ns,n_layer-1)
        r_coarse = np.random.uniform(0,6,ns+1)
        f1d = interp1d(np.linspace(0,ns,ns+1),r_coarse, kind='cubic')
        r_fine = f1d(t0)
        coef = linCoef(r_fine)
        y1 = r_fine * (1 + (np.random.random(len(r_fine)) - 0.5) * (0.015 * coef))
        x2 = np.linspace(min(t0), max(t0), 45)
        ius = InterpolatedUnivariateSpline(t0, y1)
        y2 = ius(x2)
        ius = InterpolatedUnivariateSpline(x2, y2)
        res_i = ius(t0)
    res_i[res_i < 0] = 0.0;   res_i[res_i > 6] = 6.0
    return res_i

def linCoef(T_g):
    y1 = 1
    y2 = abs(np.max(T_g)/np.min(T_g));
    coef = y1 +(y2-y1)/(np.min(T_g)-np.max(T_g))*(T_g-max(T_g))
    return coef

def rand_bin_array(K, N):
    arr = np.zeros(N)
    arr[:K]  = np.random.uniform(0,3,K)  # 1
    np.random.shuffle(arr)
    return arr

def plot_M(fn_model,res_model):
    dscale = 1.0
    m1 = Model()
    m1.read_iter_file(fn_model)
    plot_depth = m1.model_depth[1:] / dscale
    plot_model = abs(10 ** res_model[1:, 1])
    plt.semilogx(plot_model[::-1], plot_depth[::-1], ls='steps-', color='b', lw=0.8)


def smooth(x, window_len=11, window='hanning'):   # smooth test
    #  https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    """smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    see also:
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    import numpy
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s = numpy.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')
    y = numpy.convolve(w / w.sum(), s, mode='valid')   # 'same'?
    return y[(int(np.floor(window_len/2))):-int(np.floor(window_len/2))]  #check size, revised


def gauss_kern(size, sizey=None):  # smooth test
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    import scipy
    g = gauss_kern(n, sizey=ny)
    improc = scipy.signal.convolve(im,g, mode='same')
    return improc

def savgol_smooth_Jon(input0):
    from scipy.signal import savgol_filter
    iwin0 = [5,9,11,15,23]
    iorder = 3
    for iwin in iwin0:
        input0[1] = savgol_filter(input0[1].reshape(-1), iwin, iorder).reshape(-1,1)
        input0[2] = savgol_filter(input0[2].reshape(-1), iwin, iorder).reshape(-1,1)
    return input0

def savgol_smooth_Jon_0(train_data0):
    from scipy.signal import savgol_filter
    iwin0 = [5, 9, 11, 15, 23]
    iorder = 3
    if train_data0.ndim == 1:
        train_data0 = savgol_filter.reshape((1,-1))
    for i,val in enumerate(range(train_data0.shape[0])):
        for iwin in iwin0:
            train_data0[i] = savgol_filter(train_data0[i], iwin, iorder)
    return train_data0

def savgol_smooth_Jon_1(input0):
    from scipy.signal import savgol_filter
    iwin0 = [5,9,13,17,23,31,39,49,59,69]
    iorder = 3
    for iwin in iwin0:
        input0 = savgol_filter(input0, iwin, iorder)
    return input0

def savgol_smooth_Jon_2(input0):
    from scipy.signal import savgol_filter
    iwin0 = np.linspace(5,29,7,dtype=int)  #
    iorder = 3
    output = []
    for i0 in range(300):
        iwin1 = sorted(set(np.random.choice(iwin0,6)))
        temp0 = input0.reshape(-1)
        for iwin in iwin1:
            temp0 = savgol_filter(temp0, iwin, iorder)
        output.append(temp0)
    return np.array(output)

def savgol_smooth_Jon_3(input0):
    from scipy.signal import savgol_filter
    iwin0 = np.linspace(5,69,17,dtype=int)  #
    iorder = 3
    ntimes = 300
    output = []
    for i0 in range(ntimes):
        iwin1 = sorted(set(np.random.choice(iwin0,8)))
        temp0 = input0.reshape(-1)
        for iwin in iwin1:
            temp0 = savgol_filter(temp0, iwin, iorder)
        output.append(temp0)
    return np.array(output)

def phase_shift_1(phase0):
    phase0 = abs(phase0)
    for i,val in enumerate(phase0):
        if val>90:
            phase0[i] = 180-val
    return phase0

def Z_interp1d(hsi,depth,res):
    Z = np.empty(res.shape[1])
    for i in range(res.shape[1]):
        f1d = interp1d(depth,res[:,i])
        Z[i] = f1d(hsi)
    return Z

def in_hull(p, hull):
    """
    Test if points in p are within the convex hull
    """
    try:
        if not isinstance(hull, Delaunay):
            hull = Delaunay(hull)

        return hull.find_simplex(p) >= 0
    except:
        from scipy.optimize import linprog
        # Delaunay triangulation will fail if there are collinear points; in those instances
        # use linear programming (much slower) to define a convex hull.
        def in_hull_lp(points, x):
            """
            :param points:
            :param x:
            :return:
            """
            n_points = len(points)
            n_dim = len(x)
            c = np.zeros(n_points)
            A = np.r_[points.T, np.ones((1, n_points))]
            b = np.r_[x, np.ones(1)]
            lp = linprog(c, A_eq=A, b_eq=b)
            return not lp.success
        # end func
        result = []
        for cp in p:
            result.append(in_hull_lp(hull, cp))
        # end for
        return np.array(result)
