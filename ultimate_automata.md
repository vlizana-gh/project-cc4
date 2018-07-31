# Ultimate model


Discharges:
$$
w'_{ijk1} = \max(w_{ijk}-w_{(i{+}1)jk},0) V \frac{\Delta t}{\Delta x^2} + \max(+u_{ijk},0) \frac{\Delta t}{\Delta x}
$$
$$
w'_{ijk2} = \max(w_{ijk}-w_{(i{-}1)jk},0) V \frac{\Delta t}{\Delta x^2} + \max(-u_{ijk},0) \frac{\Delta t}{\Delta x}
$$
$$
w'_{ijk3} = \max(w_{ijk}-w_{i(j{+}1)k},0) V \frac{\Delta t}{\Delta x^2} + \max(+v_{ijk},0) \frac{\Delta t}{\Delta x}
$$
$$
w'_{ijk4} = \max(w_{ijk}-w_{i(j{-}1)k},0) V \frac{\Delta t}{\Delta x^2} + \max(-v_{ijk},0) \frac{\Delta t}{\Delta x}
$$

Total incoming discharge:
$$
W_{ijk} = w'_{(i{-}1)jk1} + w'_{(i{+}1)jk2} + w'_{i(j{-}1)k3} + w'_{i(j{+}1)k4}
$$

Total outgoing discharge:
$$
W'_{ijk} = w'_{ijk1} + w'_{ijk2} + w'_{ijk3} + w'_{ijk4}
$$

Momentum transference:
$$
u'_{ijkp} =  u_{ijk} \frac{d_{ijkp}}{w_{ijk}}
$$
$$
v'_{ijkp} =  v_{ijk} \frac{d_{ijkp}}{w_{ijk}}
$$

Incoming mean height change:
$$
K_{ijk1} = (w_{(i{-}1)jk} - w_{ijk}) - \frac{1}{2} W'_{(i{-}1)jk} + W'_{ijk} - \frac{1}{2} W_{ijk}
$$
$$
K_{ijk2} = (w_{(i{+}1)jk} - w_{ijk}) - \frac{1}{2} W'_{(i{+}1)jk} + W'_{ijk} - \frac{1}{2} W_{ijk}
$$
$$
K_{ijk3} = (w_{i(j{-}1)k} - w_{ijk}) - \frac{1}{2} W'_{i(j{-}1)k} + W'_{ijk} - \frac{1}{2} W_{ijk}
$$
$$
K_{ijk4} = (w_{i(j{+}1)k} - w_{ijk}) - \frac{1}{2} W'_{i(j{+}1)jk} + W'_{ijk} - \frac{1}{2} W_{ijk}
$$

Momentum generation:
$$
u''_{ijk} =
w'_{(i{-}1)jk1} \left(
    \sqrt{\max\left(0,2gK_{ijk1}+\frac{u_{(i{-}1)jk}^2}{w_{(i{-}1)jk}^2}\right)}-\frac{u_{(i{-}1)jk}}{w_{(i{-}1)jk}}
    \right)
+ w'_{(i{+}1)jk2} \left(
    -\sqrt{\max\left(0,2gK_{ijk2}+\frac{u_{(i{+}1)jk}^2}{w_{(i{+}1)jk}^2}\right)}-\frac{u_{(i{+}1)jk}}{w_{(i{+}1)jk}}
    \right)
$$

$$
v''_{ijk} =
v'_{i(j{-}1)k3} \left(
    \sqrt{\max\left(0,2gK_{ijk3}+\frac{u_{i(j{-}1)k}^2}{w_{i(j{-}1)k}^2}\right)}-\frac{u_{i(j{-}1)k}}{w_{i(j{-}1)k}}
    \right)
+ v'_{i(j{+}1)k4} \left(
    -\sqrt{\max\left(0,2gK_{ijk4}+\frac{u_{i(j{+}1)k}^2}{w_{i(j{+}1)k}^2}\right)}-\frac{u_{i(j{+}1)k}}{w_{i(j{+}1)k}}
    \right)
$$


Water update:
$$
w_{ij(k{+}1)} = w_{ijk} + w'_{(i{-}1)jk1} + w'_{(i{+}1)jk2} + w'_{i(j{-}1)k3} + w'_{i(j{+}1)k4} - \sum_{p=1}^{4} w'_{ijkp}
$$


Momentum update:
$$
u_{ij(k{+}1)} = u_{ijk} + u''_{ijk} + u'_{(i{-}1)jk1} + u'_{(i{+}1)jk2} + u'_{i(j{-}1)k3} + u'_{i(j{+}1)k4} - \sum_{p=1}^{4} u'_{ijkp}
$$

$$
v_{ij(k{+}1)} = v_{ijk} + v''_{ijk} + v'_{(i{-}1)jk1} + v'_{(i{+}1)jk2} + v'_{i(j{-}1)k3} + v'_{i(j{+}1)k4} - \sum_{p=1}^{4} v'_{ijkp}
$$
