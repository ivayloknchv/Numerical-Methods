```python
import numpy as np
```


```python
# Lagrange
```


```python
def calc_basis_poly(nodes, k, x):
    res = 1
    for i in range(len(nodes)):
        if(i!=k):
            res*=(x-nodes[i])/(nodes[k]-nodes[i])
    return res

def lagrange_poly(nodes, values, x):
    res = 0
    for i in range(len(nodes)):
        res += values[i] * calc_basis_poly(nodes, i, x)
    return res
```


```python
# Divided Differences - Normal Formula
```


```python
def divided_difference(nodes, values):
    if len(nodes) == 1:
        return values[0]
    return (divided_difference(nodes[1::], values[1::]) - divided_difference(nodes[0:-1:], values[0:-1:])) / (nodes[-1] - nodes[0])
```


```python
# Newton
```


```python
def newton_poly(nodes, values, x):
    res = 0
    multiplier = 1
    n = len(nodes)
    
    for k in range(0, n):
        current = divided_difference(nodes[0:k+1], values[0:k+1])*multiplier
        res+=current
        multiplier *= (x-nodes[k])
        
    return res   
```


```python
# Divided Difference - Extended Formula
```


```python
def divided_difference_extended(nodes, values, l, r):
    if nodes[l] == nodes[r]:
        return values[np.argmax(nodes == nodes[l]) + r - l] / math.factorial(r-l)
    return (divided_difference_extended(nodes, values, l+1, r) - divided_difference_extended(nodes, values, l, r - 1)) / (nodes[r] - nodes[l])
```


```python
# Hermite
```


```python
def hermite_poly(nodes, values, x):
    res = 0
    multiplier = 1
    
    for i in range(len(nodes)):
        res += divided_diff_extended(nodes, values, 0, i) *  multiplier
        multiplier *= (x - nodes[i])
        
    return res
```


```python
# Newton Forward TODO
```


```python
def newton_forward(n, x0, h, f):
    nodes = np.linspace(x0, x0 + n*h, n+1)
    values  = np.array(f(nodes))
```


```python
# Exponential Basis Interpolation
```


```python
def create_matrix(nodes):
    matrix = []
    
    for node in nodes:
        current_row = []

        for i in range(len(nodes)):
            current_row.append(np.e ** (i*node))
        matrix.append(current_row)

    return np.array(matrix)
```


```python
def calc_exponential_poly(poly_coeffs, x):
    res = 0
    for i in range(len(poly_coeffs)):
        res += poly_coeffs[i] * np.exp(i * x)
    return res
```


```python
# Trigonometric Basis Poly - 1st Way
```


```python
def create_matrix(nodes):
    matrix = []
    for node in nodes:
        row = [1]
        mult = 1
        for i in range(1, len(nodes)):
            if i % 2 == 1:
                row.append(np.cos(mult*node))
            else:
                row.append(np.sin(mult*node))
                mult += 1
        matrix.append(row)
    return np.array(matrix)
```


```python
def calc_trig_poly(coeffs, x):
    res = coeffs[0]
    mult = 1
    for i in range(1, len(coeffs)):
        if i % 2 == 1:
            res += coeffs[i]*np.cos(mult*x)
        else:
            res += coeffs[i]*np.sin(mult*x)
            mult += 1
    return res
```


```python
# Trigonometric Basis Poly - 2nd Way (Maybe the better one)
```


```python
def create_matrix(nodes):
    matrix = []

    for node in nodes:
        row = []
        for i in range(0,len(nodes)):
            if i == 0: 
                row.append(1)
            elif i % 2 == 1:
                row.append(np.cos((i // 2 + 1) * node))
            else:
                row.append(np.sin((i // 2) * node))
        matrix.append(row)

    return np.array(matrix)
```


```python
# Another Basis:
```


```python
def create_matrix(nodes, basis):
    matrix = []

    for node in nodes:
        row = []
        for i in range(len(nodes)):
            row.append(basis(i, node))
        matrix.append(row)
    return np.array(matrix)
```


```python
def calc_poly(coeffs, x, basis):
    res = 0
    for i, coeff in enumerate(coeffs):
        res += coeff * basis(i, x)
    return res
```


```python
# Lagrange Analogue for Trig Basis
```


```python
def calc_basis_func(nodes, k, x):
    res = 1
    for idx, node in enumerate(nodes):
        if idx != k:
            res *= np.sin((x-node)/2) / np.sin((nodes[k]-node)/2)
    return res
    
def calc_trig_poly(nodes, values, x):
    res = 0
    for i in range(len(values)):
        res = values[i] * calc_basis_func(nodes, i, x)
    return res
```


```python
# Chebyshev Nodes - 2 Ways
```


```python
n = 4 # degree of poly is 4 but we need n+1 nodes!!!!!
chebyshev_nodes = np.empty(n+1)
for k in range(n+1):
    chebyshev_nodes[k] = np.cos((2 * (k + 1) - 1) / (2 * (n+1)) * np.pi)

nodes = np.array([np.cos(((2*k - 1) * np.pi) / (2 * (n + 1))) for k in range(1, n + 2)])
```
