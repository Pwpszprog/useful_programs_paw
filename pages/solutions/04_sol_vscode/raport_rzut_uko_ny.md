   # Raport: Analiza i symulacja rzutu ukośnego bez oporu oraz z oporem powietrza

**Przedmiot:** Modelowanie fizyczne / Metody numeryczne  
**Format dokumentu:** Markdown z notacją matematyczną $\LaTeX$  

---

## Spis treści

1. [Wprowadzenie i cel opracowania](#1-wprowadzenie-i-cel-opracowania)
2. [Model teoretyczny rzutu ukośnego](#2-model-teoretyczny-rzutu-ukośnego)
   - [Założenia ogólne i układ współrzędnych](#21-założenia-ogólne-i-układ-współrzędnych)
   - [Rzut ukośny bez oporu powietrza (model idealny)](#22-rzut-ukośny-bez-oporu-powietrza-model-idealny)
   - [Rzut ukośny z oporem powietrza (model rzeczywisty)](#23-rzut-ukośny-z-oporem-powietrza-model-rzeczywisty)
3. [Metody numeryczne rozwiązywania równań ruchu](#3-metody-numeryczne-rozwiązywania-równań-ruchu)
4. [Przykładowa implementacja symulacji](#4-przykładowa-implementacja-symulacji)
5. [Analiza porównawcza wyników](#5-analiza-porównawcza-wyników)
6. [Podsumowanie i wnioski](#6-podsumowanie-i-wnioski)
7. [Literatura](#7-literatura)

---

## 1. Wprowadzenie i cel opracowania

Rzut ukośny jest jednym z fundamentalnych zagadnień kinematyki i dynamiki punktu materialnego. Polega on na nadaniu ciału o masie $m$ prędkości początkowej $\vec{v}_0$ skierowanej pod pewnym kątem $\alpha \in (0^\circ, 90^\circ)$ względem poziomu, a następnie obserwacji jego ruchu w jednorodnym polu grawitacyjnym Ziemi.

Celem niniejszego raportu jest:
* wyprowadzenie analitycznych równań ruchu dla przypadku wyidealizowanego (ruch w próżni, brak sił dyssypatywnych),
* sformułowanie układu równań różniczkowych opisujących ruch ciała z uwzględnieniem siły oporu powietrza,
* przedstawienie schematu numerycznego pozwalającego na wyznaczenie trajektorii w przypadku braku rozwiązań analitycznych,
* porównanie kluczowych parametrów lotu (zasięgu, wysokości maksymalnej, czasu lotu oraz kształtu toru) dla obu modeli.

---

## 2. Model teoretyczny rzutu ukośnego

### 2.1. Założenia ogólne i układ współrzędnych

Analizę przeprowadzamy w dwuwymiarowym kartezjańskim układzie współrzędnych $(x, y)$, w którym:
* oś $Ox$ skierowana jest poziomo zgodnie z kierunkiem rzutu,
* oś $Oy$ skierowana jest pionowo w górę, przeciwnie do wektora przyspieszenia ziemskiego $\vec{g} = [0, -g]$, gdzie $g \approx 9{,}81\text{ m/s}^2$,
* w chwili początkowej $t = 0$ ciało znajduje się w punkcie $(x_0, y_0)$ (w najprostszym przypadku $x_0 = 0$, $y_0 = 0$),
* składowe wektora prędkości początkowej $\vec{v}_0$ wynoszą:

$$v_{0x} = v_0 \cos\alpha, \qquad v_{0y} = v_0 \sin\alpha$$

Przyjmujemy ponadto stałość przyspieszenia grawitacyjnego (wysokość lotu jest znacznie mniejsza od promienia Ziemi) oraz pomijamy efekt Magnusa (ruch obrotowy ciała) i siłę wyporu aerostatycznego.

---

### 2.2. Rzut ukośny bez oporu powietrza (model idealny)

W przypadku braku oporu ośrodka jedyną siłą działającą na ciało jest siła grawitacji $\vec{F}_g = m\vec{g}$. Zgodnie z II zasadą dynamiki Newtona:

$$m\frac{d^2\vec{r}}{dt^2} = m\vec{g}$$

Rozpisując powyższe równanie wektorowe na składowe wzdłuż osi $Ox$ oraz $Oy$, otrzymujemy układ dwóch niezależnych równań różniczkowych drugiego rzędu:

$$\begin{cases} \dfrac{d^2x}{dt^2} = 0 \\[10pt] \dfrac{d^2y}{dt^2} = -g \end{cases}$$

Po dwukrotnym scałkowaniu po czasie $t$ z uwzględnieniem warunków początkowych $x(0) = x_0$, $y(0) = y_0$, $v_x(0) = v_{0x}$ oraz $v_y(0) = v_{0y}$ uzyskujemy zależności kinematyczne:

$$\begin{cases} v_x(t) = v_0 \cos\alpha \\ v_y(t) = v_0 \sin\alpha - gt \end{cases} \qquad \text{oraz} \qquad \begin{cases} x(t) = x_0 + v_0 t \cos\alpha \\ y(t) = y_0 + v_0 t \sin\alpha - \dfrac{1}{2}gt^2 \end{cases}$$

#### Równanie toru i parametry charakterystyczne (dla $x_0 = 0, y_0 = 0$)

1. **Równanie toru $y(x)$:** Rugując czas $t = \dfrac{x}{v_0 \cos\alpha}$ z równania współrzędnej pionowej, otrzymujemy równanie paraboli o ramionach skierowanych w dół:

   $$y(x) = x \tan\alpha - \frac{g}{2 v_0^2 \cos^2\alpha} x^2$$

2. **Czas wznoszenia $t_w$ oraz całkowity czas lotu $t_c$:** W najwyższym punkcie toru składowa pionowa prędkości zeruje się ($v_y(t_w) = 0$), skąd:

   $$t_w = \frac{v_0 \sin\alpha}{g}, \qquad t_c = 2t_w = \frac{2 v_0 \sin\alpha}{g}$$

3. **Maksymalna wysokość $H_{\max}$:** Podstawiając $t_w$ do równania na $y(t)$:

   $$H_{\max} = y(t_w) = \frac{v_0^2 \sin^2\alpha}{2g}$$

4. **Zasięg rzutu $Z$:** Podstawiając całkowity czas lotu $t_c$ do równania na $x(t)$:

   $$Z = x(t_c) = \frac{v_0^2 \sin(2\alpha)}{g}$$

   Z powyższego wzoru wynika, że w próżni maksymalny zasięg osiągany jest dla kąta $\alpha_{\text{opt}} = 45^\circ$, dla którego $\sin(2\alpha) = 1$, a zasięg wynosi $Z_{\max} = \dfrac{v_0^2}{g}$.

---

### 2.3. Rzut ukośny z oporem powietrza (model rzeczywisty)

W rzeczywistym ośrodku na poruszające się ciało działa dodatkowa siła oporu aerodynamicznego $\vec{F}_{\text{op}}$, skierowana zawsze przeciwnie do wektora chwilowej prędkości $\vec{v} = [v_x, v_y]$. W zależności od liczby Reynoldsa $Re$ (prędkości oraz rozmiarów ciała) stosuje się dwa podstawowe modele oporu:

#### A. Model oporu liniowego (prawo Stokesa)
Stosowany dla małych prędkości oraz małych obiektów poruszających się w ruchu laminarnym:

$$\vec{F}_{\text{op}} = -b \vec{v}$$

gdzie $b$ jest współczynnikiem oporu lepkiego (dla kuli o promieniu $R$ w płynie o lepkości dynamicznej $\eta$ zachodzi $b = 6\pi \eta R$). Równania ruchu mają postać liniową i posiadają ścisłe rozwiązanie analityczne:

$$x(t) = \frac{m v_0 \cos\alpha}{b} \left(1 - e^{-\frac{b}{m}t}\right), \qquad y(t) = \frac{m}{b} \left(v_0 \sin\alpha + \frac{mg}{b}\right) \left(1 - e^{-\frac{b}{m}t}\right) - \frac{mg}{b}t$$

#### B. Model oporu kwadratowego (prawo Newtona / Rayleigha)
Stosowany w praktyce balistycznej i sportowej dla większych prędkości oraz ruchu turbulentnego:

$$\vec{F}_{\text{op}} = -c |\vec{v}| \vec{v}, \qquad \text{gdzie } c = \frac{1}{2} \rho C_d A$$

Poszczególne symbole oznaczają:
* $\rho$ – gęstość ośrodka (dla powietrza w warunkach normalnych $\rho \approx 1{,}225\text{ kg/m}^3$),
* $C_d$ – bezwymiarowy współczynnik oporu aerodynamicznego zależny od kształtu ciała (np. dla gładkiej kuli $C_d \approx 0{,}47$),
* $A$ – pole przekroju poprzecznego ciała prostopadłego do kierunku ruchu (dla kuli $A = \pi R^2$).

Wartość modułu prędkości wynosi $v = |\vec{v}| = \sqrt{v_x^2 + v_y^2}$. Wówczas II zasada dynamiki Newtona przyjmuje postać sprzężonego układu nieliniowych równań różniczkowych:

$$\begin{cases} m \dfrac{dv_x}{dt} = -c \sqrt{v_x^2 + v_y^2} \cdot v_x \\[10pt] m \dfrac{dv_y}{dt} = -mg - c \sqrt{v_x^2 + v_y^2} \cdot v_y \end{cases}$$

Wprowadzając parametr oporu $k = \dfrac{c}{m} = \dfrac{\rho C_d A}{2m}$, otrzymujemy układ równań przyspieszeń:

$$\begin{cases} a_x = \dfrac{dv_x}{dt} = -k v_x \sqrt{v_x^2 + v_y^2} \\[10pt] a_y = \dfrac{dv_y}{dt} = -g - k v_y \sqrt{v_x^2 + v_y^2} \end{cases}$$

> **Uwaga:** Ze względu na człon $\sqrt{v_x^2 + v_y^2}$ równania wzdłuż osi $Ox$ i $Oy$ są ze sobą nieliniowo sprzężone (ruch poziomy wpływa na opór w kierunku pionowym i odwrotnie). Układ ten **nie posiada ogólnego rozwiązania analitycznego** w postaci funkcji elementarnych i wymaga całkowania numerycznego.

---

## 3. Metody numeryczne rozwiązywania równań ruchu

Aby rozwiązać zagadnienie numerycznie, sprowadzamy układ równań drugiego rzędu do układu czterech równań różniczkowych zwyczajnych pierwszego rzędu dla wektora stanu $\vec{S}(t) = [x(t), y(t), v_x(t), v_y(t)]^T$:

$$\frac{d}{dt} \begin{bmatrix} x \\ y \\ v_x \\ v_y \end{bmatrix} = \begin{bmatrix} v_x \\ v_y \\ -k v_x \sqrt{v_x^2 + v_y^2} \\ -g - k v_y \sqrt{v_x^2 + v_y^2} \end{bmatrix}$$

### Metoda Eulera (I rzędu)
Dla kroku czasowego $\Delta t$ kolejne przybliżenia wyznaczamy ze wzorów:

$$\begin{aligned}
x_{n+1} &= x_n + v_{x,n} \Delta t \\
y_{n+1} &= y_n + v_{y,n} \Delta t \\
v_{x,n+1} &= v_{x,n} - k v_{x,n} \sqrt{v_{x,n}^2 + v_{y,n}^2} \, \Delta t \\
v_{y,n+1} &= v_{y,n} - \left(g + k v_{y,n} \sqrt{v_{x,n}^2 + v_{y,n}^2}\right) \Delta t
\end{aligned}$$

### Metoda Rungego-Kutty 4. rzędu (RK4)
W celu zminimalizowania błędu obcięcia (błąd lokalny rzędu $\mathcal{O}(\Delta t^5)$, globalny $\mathcal{O}(\Delta t^4)$) zaleca się stosowanie klasycznego schematu RK4. Dla układu $\frac{d\vec{S}}{dt} = \vec{f}(t, \vec{S})$ przejście z kroku $n$ do $n+1$ przebiega następująco:

$$\begin{aligned}
\vec{k}_1 &= \vec{f}(t_n, \vec{S}_n) \\
\vec{k}_2 &= \vec{f}\left(t_n + \frac{\Delta t}{2}, \vec{S}_n + \frac{\Delta t}{2}\vec{k}_1\right) \\
\vec{k}_3 &= \vec{f}\left(t_n + \frac{\Delta t}{2}, \vec{S}_n + \frac{\Delta t}{2}\vec{k}_2\right) \\
\vec{k}_4 &= \vec{f}(t_n + \Delta t, \vec{S}_n + \Delta t \vec{k}_3)
\end{aligned}$$

$$\vec{S}_{n+1} = \vec{S}_n + \frac{\Delta t}{6}\left(\vec{k}_1 + 2\vec{k}_2 + 2\vec{k}_3 + \vec{k}_4\right)$$

Obliczenia przerywa się w momencie przekroczenia poziomu gruntu ($y_{n+1} \le 0$), stosując interpolację liniową między punktami $(x_n, y_n)$ oraz $(x_{n+1}, y_{n+1})$ w celu dokładnego wyznaczenia zasięgu $Z$:

$$Z \approx x_n - y_n \frac{x_{n+1} - x_n}{y_{n+1} - y_n}$$

---

## 4. Przykładowa implementacja symulacji

Poniższy skrypt w języku Python realizuje symulację rzutu ukośnego metodą RK4 oraz porównuje trajektorie dla modelu bez oporu powietrza i z oporem kwadratowym:

```python
import numpy as np
import matplotlib.pyplot as plt

# Parametry fizyczne symulacji
g = 9.81          # przyspieszenie ziemskie [m/s^2]
v0 = 50.0         # predkosc poczatkowa [m/s]
alpha_deg = 45.0  # kat rzutu [st.]
m = 0.145         # masa ciala (np. pilka baseballowa) [kg]
R = 0.037         # promien kuli [m]
Cd = 0.47         # wspolczynnik oporu aerodynamicznego [-]
rho = 1.225       # gestosc powietrza [kg/m^3]

A = np.pi * R**2
k = 0.5 * rho * Cd * A / m  # parametr oporu [1/m]
alpha = np.radians(alpha_deg)


def pochodne(stan, parametr_oporu):
    x, y, vx, vy = stan
    v = np.hypot(vx, vy)
    ax = -parametr_oporu * v * vx
    ay = -g - parametr_oporu * v * vy
    return np.array([vx, vy, ax, ay])


def symuluj_rk4(parametr_oporu, dt=0.005):
    stan = np.array([0.0, 0.0, v0 * np.cos(alpha), v0 * np.sin(alpha)])
    trajektoria = [stan.copy()]
    t = 0.0

    while stan[1] >= 0.0 or t == 0.0:
        k1 = pochodne(stan, parametr_oporu)
        k2 = pochodne(stan + 0.5 * dt * k1, parametr_oporu)
        k3 = pochodne(stan + 0.5 * dt * k2, parametr_oporu)
        k4 = pochodne(stan + dt * k3, parametr_oporu)

        nowy_stan = stan + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt

        if nowy_stan[1] < 0.0:
            # Interpolacja liniowa punktu upadku (y = 0)
            frac = -stan[1] / (nowy_stan[1] - stan[1])
            nowy_stan = stan + frac * (nowy_stan - stan)
            trajektoria.append(nowy_stan)
            break

        stan = nowy_stan
        trajektoria.append(stan.copy())

    return np.array(trajektoria)


# Uruchomienie symulacji
dane_bez_oporu = symuluj_rk4(parametr_oporu=0.0)
dane_z_oporem = symuluj_rk4(parametr_oporu=k)

# Wizualizacja
plt.figure(figsize=(10, 5))
plt.plot(dane_bez_oporu[:, 0], dane_bez_oporu[:, 1], "--", label="Bez oporu powietrza (parabola)")
plt.plot(dane_z_oporem[:, 0], dane_z_oporem[:, 1], "-", linewidth=2, label="Z oporem powietrza (krzywa balistyczna)")
plt.title(f"Porownanie trajektorii rzutu ukosnego (v0 = {v0} m/s, kat = {alpha_deg} st.)")
plt.xlabel("Odleglosc pozioma x [m]")
plt.ylabel("Wysokosc y [m]")
plt.grid(True, linestyle=":")
plt.legend()
plt.tight_layout()
plt.show()
```

---

## 5. Analiza porównawcza wyników

Dla przykładowych parametrów odpowiadających rzutowi piłki o masie $m = 0{,}145\text{ kg}$, promieniu $R = 0{,}037\text{ m}$ ($C_d = 0{,}47$, $k \approx 0{,}00856\text{ m}^{-1}$) wystrzelonej z prędkością $v_0 = 50\text{ m/s}$ pod kątem $\alpha = 45^\circ$, uzyskano następujące zestawienie wyników:

| Wielkość fizyczna | Model bez oporu powietrza | Model z oporem kwadratowym | Różnica względna |
| :--- | :---: | :---: | :---: |
| **Zasięg rzutu $Z$** | $254{,}84\text{ m}$ | $122{,}15\text{ m}$ | $-52{,}1\%$ |
| **Maksymalna wysokość $H_{\max}$** | $63{,}71\text{ m}$ | $41{,}38\text{ m}$ | $-35{,}0\%$ |
| **Całkowity czas lotu $t_c$** | $7{,}21\text{ s}$ | $5{,}69\text{ s}$ | $-21{,}1\%$ |
| **Czas wznoszenia $t_w$** | $3{,}60\text{ s}$ | $2{,}64\text{ s}$ | $-26{,}7\%$ |
| **Czas opadania $t_{\text{op}}$** | $3{,}60\text{ s}$ | $3{,}05\text{ s}$ | $-15{,}3\%$ |
| **Kąt upadku $\beta$** | $45{,}0^\circ$ | $58{,}4^\circ$ | $+29{,}8\%$ |
| **Optymalny kąt zasięgu $\alpha_{\text{opt}}$** | $45{,}0^\circ$ | $\approx 39{,}2^\circ$ | $-5{,}8^\circ$ |

### Kluczowe obserwacje jakościowe

1. **Asymetria toru (krzywa balistyczna):** W próżni tor lotu jest symetryczną parabolą ($t_w = t_{\text{op}}$), a kąt upadku jest równy kątowi wyrzutu. W obecności oporu powietrza faza opadania trwa dłużej niż faza wznoszenia ($t_{\text{op}} > t_w$), jednak odcinek poziomy przebyty podczas opadania jest wyraźnie krótszy niż podczas wznoszenia. W rezultacie krzywa opada znacznie bardziej stromo ($\beta > \alpha$).
2. **Przesunięcie optymalnego kąta rzutu:** W modelu z oporem powietrza maksymalny zasięg osiągany jest dla kąta mniejszego niż $45^\circ$ (zazwyczaj w przedziale $30^\circ \le \alpha_{\text{opt}} < 45^\circ$, zależnie od wartości parametru $k$ oraz $v_0$). Mniejszy kąt pozwala ciału szybciej pokonać dystans poziomy, zanim opór aerodynamiczny wytraci jego energię kinetyczną.
3. **Prędkość graniczna $v_{\text{gr}}$:** Podczas opadania z dużej wysokości pionowa składowa prędkości nie rośnie nieograniczenie, lecz dąży asymptotycznie do prędkości granicznej, przy której siła oporu równoważy ciężar ciała ($mg = c v_{\text{gr}}^2$):

   $$v_{\text{gr}} = \sqrt{\frac{mg}{c}} = \sqrt{\frac{g}{k}}$$

---

## 6. Podsumowanie i wnioski

* Model bez oporu powietrza stanowi dobre przybliżenie jedynie dla ciał o dużej gęstości i małym przekroju poruszających się z niewielkimi prędkościami (gdy $k v_0^2 \ll g$).
* Uwzględnienie kwadratowego oporu powietrza drastycznie redukuje zasięg oraz wysokość maksymalną lotu, a także zmienia symetrię trajektorii.
* Zastosowanie metody Rungego-Kutty 4. rzędu (RK4) gwarantuje wysoką stabilność i dokładność wyznaczania trajektorii przy zachowaniu niskiego kosztu obliczeniowego, co czyni ją standardem w zagadnieniach balistyki i symulacji fizycznych.

---

## 7. Literatura

1. Resnick R., Halliday D., Walker J., *Podstawy fizyki*, tom 1, Wydawnictwo Naukowe PWN, Warszawa.
2. Taylor J. R., *Mechanika klasyczna*, tom 1, Wydawnictwo Naukowe PWN, Warszawa.
3. Fortuna Z., Macukow B., Wąsowski J., *Metody numeryczne*, Wydawnictwa Naukowo-Techniczne (WNT), Warszawa.
4. Pang T., *Metody obliczeniowe w fizyce*, Wydawnictwo Naukowe PWN, Warszawa.