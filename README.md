# CC4 - Flooding Model

## Phase 1.- Cellular Automata

The water transference between two adjacent cells $i$ and $j$ is:
$$ Q_{ij} = 
\begin{cases}
-\alpha \min(\omega_i, (t_i + w_i) - (t_j + w_j)) \quad & \if t_i + \omega_i \geqslant t_j + \omega_j\\
+\alpha \min(\omega_j, (t_j + w_j) - (t_i + w_i)) \quad & \if t_i + \omega_i \geqslant t_j + \omega_j
\end{cases} $$

where $\alpha$ is a constant in $[0,1]$.
