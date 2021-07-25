+++
title = "Multivariate Bernoulli Distribution: Bayesian Analysis with PyMC3"
date = 2021-07-25T20:15:00+05:30
draft = false

# Tags and categories
# For example, use `tags = []` for no tags, or the form `tags = ["A Tag", "Another Tag"]` for one or more tags.
tags = ["bernoulli", "multivariate bernoulli", "bayesian", "pymc3", "python"]
categories = ["Statistics", "Bayesian"]

# Featured image
# Place your image in the `static/img/` folder and reference its filename below, e.g. `image = "example.jpg"`.
# Use `caption` to display an image caption.
#   Markdown linking is allowed, e.g. `caption = "[Image credit](http://example.org)"`.
# Set `preview` to `false` to disable the thumbnail in listings.
[header]
image = ""
caption = ""
preview = true

+++



```python
import matplotlib.pyplot as plt
import pymc3 as pm
import aesara.tensor as at
import arviz as az
import scipy.stats as st
import numpy as np
```


```python
az.style.use('arviz-whitegrid')
```

## Multivariate bernoulli distribution

### Sources

1. https://doi.org/10.1016/0047-259X(90)90084-U

2. https://arxiv.org/pdf/1206.1874

### Motivation

I have recently been introduced to Bayesian framework for statistical analysis with [PyMC3 package](https://docs.pymc.io/) for Probabilistic Programming in Python 3.

I came across to a data given to me by my friend, which had multiple columns of binary data (responses to various questions with answers in binary format, yes or no). Each column could be represented by bernoulli process with an unknown probability of getting answer "yes", which is what I wished to estimate. In addition to the above, the data also contains possible association between the columns, which I so dearly wanted to estimate. After a lot of literature search I came across the work done by B Dai in 2012 - 2013, which have been cited above in Sources.

The following post aims to translate the log likelihood function, joint probability distribution, marginal probability distribution and conditions for independence described above into PyMC3 code, which can be used by others to solve related problems.

### Formulation

Let \\(Y = (Y\_1, Y\_2, ..., Y\_K)\\) be a \\(K\\) dimensional random vector of possibly correlated Bernoulli random variables and let \\(y = (y\_1, y\_2, ..., y\_K)\\) be a realisation of \\(Y\\). The more general form of \\(p(y\_1, y\_2, ..., y\_K)\\) of joint probability density is

\\[
    P(Y\_1 = y\_1, ..., Y\_K = y\_K) = p\_{0,0,..,0}^{\\prod\_{j=1}^{K}{(1-y\_j)}} p\_{1,0,..,0}^{y\_1\\prod\_{j=2}^{K}{(1-y\_j)}} p\_{0,1,..,0}^{y\_1y\_2\\prod\_{j=3}^{K}{(1-y\_j)}} p\_{1,1,..,1}^{\\prod\_{j=}^{K}{y\_j}}
\\]

To simplify the notation, denote quantity \\(S\\) to be

\\[
S^{j\_1j\_2...j\_r} = \\sum\_{1 \\leq s \\leq r}{f^{j\_s}} + \\sum\_{1 \\leq s < t \\leq r}{f^{j\_sj\_t}} + ... + f^{j\_1j\_2...j\_r}
\\]

Also, we will define the interaction function \\(B\\) as

\\[
B^{j\_1j\_2...j\_r}(y) = y\_{j\_1}y\_{j\_2}...y\_{j\_r}
\\]

The log-linear formulation of multivariate Bernoulli distribution will be

\\[
l(y, \\mathbf{f}) = -log[p(y)] = -\\left[\\sum\_{r=1}^{K}{\\left(\\sum\_{1 \\leq j\_1 < j\_2 < ... < j\_r \\leq K}{f^{j\_1j\_2...j\_r}B^{j\_1j\_2...j\_r}(y)} \\right)} - b(\\mathbf{f})\\right]
\\]

where, \\(\\mathbf{f} = (f^1, f^2, ..., f^{12...K})^T\\), is the vector of natural parameters for the multivariate Bernoulli distribution.

The normalising factor \\(b(\\mathbf{f}) = p\_{0,0,...,0}\\) is defined as

\\[
b(\\mathbf{f}) = log \\sum\_{r=1}^{K}{\\left [1 + \\left (\\sum\_{1 \\leq j\_1 < j\_2 < ... < j\_r \\leq K}{e^{S^{j\_1j\_2...j\_r}}}\\right ) \\right ]}
\\]

### PyMC3 implementation


```python
from itertools import chain, combinations

def powerset(iterable):
    '''
    [0,1,2] -> [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
    '''
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))
```


```python
def bernoulli_mv(f, b_f, subsets):
    '''
    Multivariate Bernoulli Distribution
    ===================================
    f: vector of f
    subsets: list of tuples of indices like [0,1,2] -> [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
    '''
    
    def logp_(y):
        res = 0
        for i, sub in enumerate(subsets):
            res += pm.math.sum(f[i] * pm.math.prod(y[:, sub], axis = 1))
        
        res = res - (b_f * y.shape[0])
        
        return res
    
    return logp_
```


```python
y = np.random.choice([0, 1], (50,5))
```


```python
subsets = list((powerset(range(y.shape[1]))))
```


```python
subsets
```




    [(0,),
     (1,),
     (2,),
     (3,),
     (4,),
     (0, 1),
     (0, 2),
     (0, 3),
     (0, 4),
     (1, 2),
     (1, 3),
     (1, 4),
     (2, 3),
     (2, 4),
     (3, 4),
     (0, 1, 2),
     (0, 1, 3),
     (0, 1, 4),
     (0, 2, 3),
     (0, 2, 4),
     (0, 3, 4),
     (1, 2, 3),
     (1, 2, 4),
     (1, 3, 4),
     (2, 3, 4),
     (0, 1, 2, 3),
     (0, 1, 2, 4),
     (0, 1, 3, 4),
     (0, 2, 3, 4),
     (1, 2, 3, 4),
     (0, 1, 2, 3, 4)]




```python
with pm.Model() as mod:
    fs = pm.Normal('f', 0, 2, shape = len(subsets))

    b_f = 1
    for s in subsets:
        ss = list(powerset(s))
        b_f += pm.math.exp(pm.math.sum(fs[np.array([s in ss for s in subsets])]))
    b_f = pm.Deterministic('b_f', pm.math.log(b_f))
    
    ys = pm.DensityDist('ys', bernoulli_mv(fs, b_f, subsets), observed = {'y': y})
```


```python
with mod:
    trace = pm.sample(return_inferencedata=True)
```

    Auto-assigning NUTS sampler...
    Initializing NUTS using jitter+adapt_diag...
    /home/suman/miniconda3/lib/python3.8/site-packages/aesara/graph/fg.py:525: UserWarning: Variable Elemwise{mul,no_inplace}.0 cannot be replaced; it isn't in the FunctionGraph
      warnings.warn(
    Multiprocess sampling (2 chains in 2 jobs)
    NUTS: [f]




<div>
    <style>
        /* Turns off some styling */
        progress {
            /* gets rid of default border in Firefox and Opera. */
            border: none;
            /* Needs to be in here for Safari polyfill so background images work as expected. */
            background-size: auto;
        }
        .progress-bar-interrupted, .progress-bar-interrupted::-webkit-progress-bar {
            background: #F44336;
        }
    </style>
  <progress value='4000' class='' max='4000' style='width:300px; height:20px; vertical-align: middle;'></progress>
  100.00% [4000/4000 02:41<00:00 Sampling 2 chains, 0 divergences]
</div>



    Sampling 2 chains for 1_000 tune and 1_000 draw iterations (2_000 + 2_000 draws total) took 166 seconds.



```python
az.plot_trace(trace)
```

```python
sumry = az.summary(trace)
sumry.index = subsets + ['b_f']
sumry
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>mean</th>
      <th>sd</th>
      <th>hdi_3%</th>
      <th>hdi_97%</th>
      <th>mcse_mean</th>
      <th>mcse_sd</th>
      <th>ess_bulk</th>
      <th>ess_tail</th>
      <th>r_hat</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>(0,)</th>
      <td>0.105</td>
      <td>0.811</td>
      <td>-1.492</td>
      <td>1.486</td>
      <td>0.026</td>
      <td>0.018</td>
      <td>982.0</td>
      <td>1376.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(1,)</th>
      <td>-0.224</td>
      <td>0.807</td>
      <td>-1.653</td>
      <td>1.322</td>
      <td>0.030</td>
      <td>0.021</td>
      <td>705.0</td>
      <td>796.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(2,)</th>
      <td>-1.284</td>
      <td>0.867</td>
      <td>-3.030</td>
      <td>0.183</td>
      <td>0.031</td>
      <td>0.022</td>
      <td>779.0</td>
      <td>995.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(3,)</th>
      <td>0.205</td>
      <td>0.737</td>
      <td>-1.250</td>
      <td>1.529</td>
      <td>0.025</td>
      <td>0.017</td>
      <td>891.0</td>
      <td>1088.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(4,)</th>
      <td>1.300</td>
      <td>0.703</td>
      <td>-0.049</td>
      <td>2.620</td>
      <td>0.030</td>
      <td>0.021</td>
      <td>538.0</td>
      <td>1077.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1)</th>
      <td>0.846</td>
      <td>0.945</td>
      <td>-0.828</td>
      <td>2.673</td>
      <td>0.031</td>
      <td>0.023</td>
      <td>917.0</td>
      <td>1140.0</td>
      <td>1.01</td>
    </tr>
    <tr>
      <th>(0, 2)</th>
      <td>0.678</td>
      <td>1.003</td>
      <td>-1.277</td>
      <td>2.468</td>
      <td>0.034</td>
      <td>0.024</td>
      <td>887.0</td>
      <td>927.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 3)</th>
      <td>0.012</td>
      <td>0.949</td>
      <td>-1.762</td>
      <td>1.831</td>
      <td>0.033</td>
      <td>0.023</td>
      <td>833.0</td>
      <td>1173.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 4)</th>
      <td>-2.639</td>
      <td>1.109</td>
      <td>-4.892</td>
      <td>-0.755</td>
      <td>0.034</td>
      <td>0.024</td>
      <td>1075.0</td>
      <td>1355.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(1, 2)</th>
      <td>0.588</td>
      <td>0.975</td>
      <td>-1.185</td>
      <td>2.402</td>
      <td>0.031</td>
      <td>0.023</td>
      <td>978.0</td>
      <td>1196.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(1, 3)</th>
      <td>-0.210</td>
      <td>0.962</td>
      <td>-2.129</td>
      <td>1.499</td>
      <td>0.039</td>
      <td>0.028</td>
      <td>610.0</td>
      <td>870.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(1, 4)</th>
      <td>0.072</td>
      <td>0.908</td>
      <td>-1.601</td>
      <td>1.734</td>
      <td>0.034</td>
      <td>0.024</td>
      <td>721.0</td>
      <td>1003.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(2, 3)</th>
      <td>0.721</td>
      <td>0.990</td>
      <td>-1.192</td>
      <td>2.561</td>
      <td>0.033</td>
      <td>0.026</td>
      <td>878.0</td>
      <td>1012.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(2, 4)</th>
      <td>0.657</td>
      <td>0.960</td>
      <td>-1.319</td>
      <td>2.210</td>
      <td>0.034</td>
      <td>0.024</td>
      <td>818.0</td>
      <td>1054.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(3, 4)</th>
      <td>-0.691</td>
      <td>0.865</td>
      <td>-2.355</td>
      <td>0.923</td>
      <td>0.027</td>
      <td>0.019</td>
      <td>1031.0</td>
      <td>947.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1, 2)</th>
      <td>-0.950</td>
      <td>1.164</td>
      <td>-3.118</td>
      <td>1.183</td>
      <td>0.035</td>
      <td>0.026</td>
      <td>1087.0</td>
      <td>920.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1, 3)</th>
      <td>-0.967</td>
      <td>1.189</td>
      <td>-3.129</td>
      <td>1.154</td>
      <td>0.040</td>
      <td>0.028</td>
      <td>878.0</td>
      <td>948.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1, 4)</th>
      <td>0.432</td>
      <td>1.195</td>
      <td>-1.754</td>
      <td>2.582</td>
      <td>0.036</td>
      <td>0.026</td>
      <td>1084.0</td>
      <td>1280.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 2, 3)</th>
      <td>-0.200</td>
      <td>1.171</td>
      <td>-2.328</td>
      <td>1.996</td>
      <td>0.037</td>
      <td>0.027</td>
      <td>974.0</td>
      <td>1180.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 2, 4)</th>
      <td>0.579</td>
      <td>1.259</td>
      <td>-1.800</td>
      <td>2.954</td>
      <td>0.040</td>
      <td>0.028</td>
      <td>995.0</td>
      <td>1087.0</td>
      <td>1.01</td>
    </tr>
    <tr>
      <th>(0, 3, 4)</th>
      <td>-0.727</td>
      <td>1.414</td>
      <td>-3.555</td>
      <td>1.626</td>
      <td>0.042</td>
      <td>0.036</td>
      <td>1139.0</td>
      <td>943.0</td>
      <td>1.01</td>
    </tr>
    <tr>
      <th>(1, 2, 3)</th>
      <td>1.388</td>
      <td>1.144</td>
      <td>-0.822</td>
      <td>3.425</td>
      <td>0.038</td>
      <td>0.027</td>
      <td>927.0</td>
      <td>1070.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(1, 2, 4)</th>
      <td>-0.323</td>
      <td>1.125</td>
      <td>-2.514</td>
      <td>1.715</td>
      <td>0.034</td>
      <td>0.024</td>
      <td>1076.0</td>
      <td>1242.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(1, 3, 4)</th>
      <td>0.622</td>
      <td>1.099</td>
      <td>-1.404</td>
      <td>2.655</td>
      <td>0.037</td>
      <td>0.026</td>
      <td>861.0</td>
      <td>1056.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(2, 3, 4)</th>
      <td>-0.930</td>
      <td>1.132</td>
      <td>-3.074</td>
      <td>1.220</td>
      <td>0.034</td>
      <td>0.024</td>
      <td>1129.0</td>
      <td>1156.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1, 2, 3)</th>
      <td>-0.570</td>
      <td>1.340</td>
      <td>-3.090</td>
      <td>1.855</td>
      <td>0.040</td>
      <td>0.028</td>
      <td>1132.0</td>
      <td>1169.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1, 2, 4)</th>
      <td>1.173</td>
      <td>1.407</td>
      <td>-1.447</td>
      <td>3.790</td>
      <td>0.043</td>
      <td>0.032</td>
      <td>1080.0</td>
      <td>1207.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1, 3, 4)</th>
      <td>1.452</td>
      <td>1.487</td>
      <td>-1.266</td>
      <td>4.308</td>
      <td>0.047</td>
      <td>0.038</td>
      <td>1016.0</td>
      <td>1039.0</td>
      <td>1.01</td>
    </tr>
    <tr>
      <th>(0, 2, 3, 4)</th>
      <td>0.101</td>
      <td>1.554</td>
      <td>-2.724</td>
      <td>3.120</td>
      <td>0.045</td>
      <td>0.040</td>
      <td>1219.0</td>
      <td>1013.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(1, 2, 3, 4)</th>
      <td>-2.363</td>
      <td>1.360</td>
      <td>-4.810</td>
      <td>0.373</td>
      <td>0.043</td>
      <td>0.032</td>
      <td>1018.0</td>
      <td>1085.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>(0, 1, 2, 3, 4)</th>
      <td>1.404</td>
      <td>1.547</td>
      <td>-1.584</td>
      <td>4.178</td>
      <td>0.039</td>
      <td>0.034</td>
      <td>1544.0</td>
      <td>1164.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>b_f</th>
      <td>3.955</td>
      <td>0.622</td>
      <td>2.760</td>
      <td>5.096</td>
      <td>0.022</td>
      <td>0.015</td>
      <td>832.0</td>
      <td>982.0</td>
      <td>1.00</td>
    </tr>
  </tbody>
</table>
</div>



### f \\(\\rightarrow\\) Joint Probabilities

\\[
p(j\_1, j\_2, ..., j\_r\\ positions\\ are\\ 1,\\ others\\ are\\ 0) = \\frac{e^{S^{j\_1j\_2...j\_r}}}{e^{b(\\mathbf{f})}}
\\]

We can get all the joint probabilities as follows


```python
p_names = []
for sub in subsets:
    p_name = ['0'] * y.shape[1]
    for i in sub:
        p_name[i] = '1'
    p_names.append(p_name)

p_names = ['p_' + ''.join(name) for name in p_names]

for p_name, s in zip(p_names, subsets):
    ss = list(powerset(s))
    idx_S = np.array([s in ss for s in subsets])
    trace.posterior[p_name] = np.exp(np.sum(trace.posterior['f'][...,idx_S], axis = 2)) / \
                              np.exp(trace.posterior['b_f'])

trace.posterior['p_' + '0' * y.shape[1]] = np.exp(-trace.posterior['b_f'])
```


```python
az.summary(trace, var_names = ['p_'], filter_vars='like')
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>mean</th>
      <th>sd</th>
      <th>hdi_3%</th>
      <th>hdi_97%</th>
      <th>mcse_mean</th>
      <th>mcse_sd</th>
      <th>ess_bulk</th>
      <th>ess_tail</th>
      <th>r_hat</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>p_10000</th>
      <td>0.027</td>
      <td>0.018</td>
      <td>0.003</td>
      <td>0.063</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1269.0</td>
      <td>1114.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01000</th>
      <td>0.019</td>
      <td>0.014</td>
      <td>0.002</td>
      <td>0.044</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>932.0</td>
      <td>1178.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00100</th>
      <td>0.008</td>
      <td>0.007</td>
      <td>0.000</td>
      <td>0.020</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>780.0</td>
      <td>1120.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00010</th>
      <td>0.028</td>
      <td>0.018</td>
      <td>0.005</td>
      <td>0.060</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1356.0</td>
      <td>1359.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00001</th>
      <td>0.077</td>
      <td>0.034</td>
      <td>0.022</td>
      <td>0.142</td>
      <td>0.001</td>
      <td>0.001</td>
      <td>1559.0</td>
      <td>1458.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11000</th>
      <td>0.047</td>
      <td>0.027</td>
      <td>0.006</td>
      <td>0.095</td>
      <td>0.001</td>
      <td>0.000</td>
      <td>1579.0</td>
      <td>1282.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_10100</th>
      <td>0.016</td>
      <td>0.014</td>
      <td>0.001</td>
      <td>0.043</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1319.0</td>
      <td>1217.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_10010</th>
      <td>0.033</td>
      <td>0.021</td>
      <td>0.002</td>
      <td>0.068</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1837.0</td>
      <td>1401.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_10001</th>
      <td>0.009</td>
      <td>0.010</td>
      <td>0.000</td>
      <td>0.026</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1430.0</td>
      <td>1357.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01100</th>
      <td>0.012</td>
      <td>0.011</td>
      <td>0.000</td>
      <td>0.030</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1416.0</td>
      <td>1003.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01010</th>
      <td>0.021</td>
      <td>0.016</td>
      <td>0.000</td>
      <td>0.050</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1312.0</td>
      <td>1325.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01001</th>
      <td>0.068</td>
      <td>0.032</td>
      <td>0.013</td>
      <td>0.124</td>
      <td>0.001</td>
      <td>0.001</td>
      <td>1914.0</td>
      <td>1092.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00110</th>
      <td>0.019</td>
      <td>0.015</td>
      <td>0.001</td>
      <td>0.046</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1615.0</td>
      <td>1232.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00101</th>
      <td>0.045</td>
      <td>0.027</td>
      <td>0.006</td>
      <td>0.092</td>
      <td>0.001</td>
      <td>0.000</td>
      <td>1930.0</td>
      <td>1295.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00011</th>
      <td>0.051</td>
      <td>0.029</td>
      <td>0.006</td>
      <td>0.104</td>
      <td>0.001</td>
      <td>0.000</td>
      <td>1856.0</td>
      <td>1548.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11100</th>
      <td>0.022</td>
      <td>0.018</td>
      <td>0.001</td>
      <td>0.055</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1921.0</td>
      <td>1436.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11010</th>
      <td>0.021</td>
      <td>0.018</td>
      <td>0.000</td>
      <td>0.052</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1718.0</td>
      <td>1417.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11001</th>
      <td>0.024</td>
      <td>0.019</td>
      <td>0.001</td>
      <td>0.057</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1884.0</td>
      <td>1303.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_10110</th>
      <td>0.031</td>
      <td>0.022</td>
      <td>0.001</td>
      <td>0.070</td>
      <td>0.001</td>
      <td>0.000</td>
      <td>1735.0</td>
      <td>1441.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_10101</th>
      <td>0.016</td>
      <td>0.015</td>
      <td>0.000</td>
      <td>0.044</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1869.0</td>
      <td>1480.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_10011</th>
      <td>0.004</td>
      <td>0.007</td>
      <td>0.000</td>
      <td>0.014</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1274.0</td>
      <td>913.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01110</th>
      <td>0.071</td>
      <td>0.035</td>
      <td>0.011</td>
      <td>0.130</td>
      <td>0.001</td>
      <td>0.001</td>
      <td>2181.0</td>
      <td>1661.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01101</th>
      <td>0.050</td>
      <td>0.028</td>
      <td>0.006</td>
      <td>0.100</td>
      <td>0.001</td>
      <td>0.000</td>
      <td>2134.0</td>
      <td>1634.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01011</th>
      <td>0.064</td>
      <td>0.033</td>
      <td>0.014</td>
      <td>0.125</td>
      <td>0.001</td>
      <td>0.001</td>
      <td>2077.0</td>
      <td>1318.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00111</th>
      <td>0.026</td>
      <td>0.020</td>
      <td>0.001</td>
      <td>0.063</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>2007.0</td>
      <td>1581.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11110</th>
      <td>0.030</td>
      <td>0.022</td>
      <td>0.002</td>
      <td>0.072</td>
      <td>0.001</td>
      <td>0.000</td>
      <td>1971.0</td>
      <td>1584.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11101</th>
      <td>0.061</td>
      <td>0.033</td>
      <td>0.010</td>
      <td>0.121</td>
      <td>0.001</td>
      <td>0.001</td>
      <td>1785.0</td>
      <td>1539.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11011</th>
      <td>0.020</td>
      <td>0.018</td>
      <td>0.000</td>
      <td>0.053</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1957.0</td>
      <td>1526.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_10111</th>
      <td>0.006</td>
      <td>0.009</td>
      <td>0.000</td>
      <td>0.022</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>1804.0</td>
      <td>1358.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_01111</th>
      <td>0.018</td>
      <td>0.017</td>
      <td>0.000</td>
      <td>0.050</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>2120.0</td>
      <td>1667.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_11111</th>
      <td>0.033</td>
      <td>0.023</td>
      <td>0.002</td>
      <td>0.077</td>
      <td>0.001</td>
      <td>0.000</td>
      <td>1862.0</td>
      <td>1732.0</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>p_00000</th>
      <td>0.023</td>
      <td>0.014</td>
      <td>0.003</td>
      <td>0.049</td>
      <td>0.000</td>
      <td>0.000</td>
      <td>832.0</td>
      <td>982.0</td>
      <td>1.01</td>
    </tr>
  </tbody>
</table>
</div>



### Marginal Probabilities

We represent marginal probability as follows

\\(p\_i\\) is probability of \\(i^{th}\\) column being 1, where \\(i \\in \\{1,2,...,K\\}\\)

Similarly, \\(p\_{ij}\\) is the probability of \\(i^{th}\\) and \\(j^{th}\\) columns being 1, where \\(i, j \\in \\{1,2,...,K\\}\\).

\\[
p\_i = \\sum\_{for\\ all\\ j\_i = 1}{p\_{j\_1j\_2...j\_k}}
\\]


```python
def marg_prob(idxs, trace, num_columns):
    '''
    idxs: list of indexes which are 1. eg [0,1] means first and second columns are 1
    trace: trace of pymc3 model output
    '''
    assert max(idxs) < num_columns
    keys = [i for i in trace.posterior.keys() if i.startswith('p_')]
    res = []
    for k in keys:
        k_num = k[2:]
        for c, i in enumerate(idxs):
            if k_num[i] != '1':
                break
            if c + 1 == len(idxs):
                res.append('p_' + k_num)
    return np.sum(trace.posterior[res].to_array(), axis = 0)
```


```python
for i in range(y.shape[1]):
    print(f'Marginal Probability {i + 1} ---')
    print(az.summary(marg_prob([i], trace, y.shape[1]), kind = 'stats'))
```

    Marginal Probability 1 ---
        mean     sd  hdi_3%  hdi_97%
    x  0.573  0.067   0.442     0.69
    Marginal Probability 2 ---
       mean     sd  hdi_3%  hdi_97%
    x   0.4  0.065    0.27    0.511
    Marginal Probability 3 ---
       mean     sd  hdi_3%  hdi_97%
    x  0.58  0.067    0.46    0.706
    Marginal Probability 4 ---
        mean     sd  hdi_3%  hdi_97%
    x  0.463  0.071   0.323    0.588
    Marginal Probability 5 ---
        mean    sd  hdi_3%  hdi_97%
    x  0.477  0.07   0.345    0.603


### Independence of outcomes

The outcomes are independent, if

\\[
S^{j\_1j\_2...j\_r} - \\sum\_{k = 1}^{r}{f^{j\_k}} = 0,\\ \\forall\\ r \\geq 2
\\]


```python
def trace_of_independence(idxs, trace, subsets, num_columns):
    '''
        idxs: list of indexes whose independence is to be checked, eg [0,1,2] means that 
        independence of 1st, 2nd and 3rd index needs to be checked
        trace: pymc3 trace, inferenceData
        subsets: list depicting powerset of the indexes
        num_columns: number of outcomes
    '''
    assert len(idxs) > 1
    
    assert max(idxs) < num_columns
    
    ss = [i for i in powerset(idxs) if len(i) > 1]
    idx_S = np.array([s in ss for s in subsets])
    
    return np.sum(trace.posterior['f'][...,idx_S], axis = 2)
```

### Issues

The sampling takes prohibitively long time when the number of variables \\(K\\) exceeds 7 - 8. I am in the process to tackle the issue.

Also, I have to make mechanism by which random samples can be generated.
