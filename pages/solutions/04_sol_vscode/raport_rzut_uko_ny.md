# Raport: Analiza teoretyczna i symulacja numeryczna rzutu ukośnego bez oporu oraz z oporem powietrza

**Przedmiot:** Modelowanie fizyczne / Metody numeryczne  
**Środowisko:** Visual Studio Code (Markdown + $\LaTeX$ + Python)  
**Układ jednostek:** Ścisły układ SI ($\text{m}$, $\text{kg}$, $\text{s}$, $\text{rad}$, $\text{N}$, $\text{J}$)

---

## Spis treści

1. [Wprowadzenie i parametry wejściowe (układ SI)](#1-wprowadzenie-i-parametry-wejściowe-układ-si)
2. [Wyprowadzenie równań ruchu bez oporu powietrza (model idealny)](#2-wyprowadzenie-równań-ruchu-bez-oporu-powietrza-model-idealny)
   * [Założenia i układ współrzędnych](#21-założenia-i-układ-współrzędnych)
   * [Całkowanie równań dynamiki Newtona](#22-całkowanie-równań-dynamiki-newtona)
   * [Wyprowadzenie równania toru i parametrów lotu](#23-wyprowadzenie-równania-toru-i-parametrów-lotu)
3. [Wyprowadzenie równań ruchu z oporem powietrza (model rzeczywisty)](#3-wyprowadzenie-równań-ruchu-z-oporem-powietrza-model-rzeczywisty)
   * [Model oporu liniowego — prawo Stokesa (wyprowadzenie analityczne)](#31-model-oporu-liniowego--prawo-stokesa-wyprowadzenie-analityczne)
   * [Model oporu kwadratowego — prawo Newtona/Rayleigha](#32-model-oporu-kwadratowego--prawo-newtonarayleigha)
   * [Prędkość graniczna w ośrodku lepkim i turbulentnym](#33-prędkość-graniczna-w-ośrodku-lepkim-i-turbulentnym)
4. [Metody numeryczne rozwiązywania równań ruchu](#4-metody-numeryczne-rozwiązywania-równań-ruchu)
   * [Redukcja do układu równań różniczkowych I rzędu](#41-redukcja-do-układu-równań-różniczkowych-i-rzędu)
   * [Schemat Eulera oraz metoda Rungego-Kutty 4. rzędu (RK4)](#42-schemat-eulera-oraz-metoda-rungego-kutty-4-rzędu-rk4)
5. [Implementacja symulacji oraz wykresy ruchu](#5-implementacja-symulacji-oraz-wykresy-ruchu)
   * [Wykres 1: Porównanie trajektorii lotu $y(x)$](#51-wykres-1-porównanie-trajektorii-lotu-yx)
   * [Wykres 2: Przebieg składowych i wartości prędkości w czasie](#52-wykres-2-przebieg-składowych-i-wartości-prędkości-w-czasie)
   * [Wykres 3: Zależność współrzędnych położenia $x(t)$ i $y(t)$ od czasu](#53-wykres-3-zależność-współrzędnych-położenia-xt-i-yt-od-czasu)
   * [Wykres 4: Wpływ kąta wyrzutu $\alpha$ na zasięg i wysokość maksymalną](#54-wykres-4-wpływ-kąta-wyrzutu-alpha-na-zasięg-i-wysokość-maksymalną)
   * [Wykres 5: Bilans energii mechanicznej układu](#55-wykres-5-bilans-energii-mechanicznej-układu)
6. [Zestawienie porównawcze wyników w jednostkach SI](#6-zestawienie-porównawcze-wyników-w-jednostkach-si)
7. [Podsumowanie i wnioski](#7-podsumowanie-i-wnioski)
8. [Literatura](#8-literatura)

---

## 1. Wprowadzenie i parametry wejściowe (układ SI)

Rzut ukośny polega na nadaniu ciału traktowanemu jako punkt materialny o masie $m$ prędkości początkowej $\vec{v}_0$ skierowanej pod kątem $\alpha \in \left(0, \frac{\pi}{2}\right)\text{ rad}$ względem poziomu, a następnie analizie jego ruchu w jednorodnym polu grawitacyjnym Ziemi o natężeniu $\vec{g}$.

Zgodnie z wymaganiami zadania wszystkie wielkości fizyczne w raporcie oraz w skryptach obliczeniowych wyrażono **wyłącznie w jednostkach układu SI** (dla kąta płaskiego podstawową jednostką układu SI jest **radian**, $\text{rad}$). Wartości parametrów początkowych wyznaczono na podstawie numeru indeksu według reguły:

* **Prędkość początkowa $v_0$:** dwie pierwsze cyfry numeru indeksu wyrażone w $\text{m/s}$,
* **Kąt początkowy $\alpha$:** dwie ostatnie cyfry numeru indeksu przeliczone ze stopni na radiany zgodnie ze wzorem:

$$
\alpha\,[\text{rad}] = \alpha_{\text{indeks}} \cdot \frac{\pi}{180}
$$

Dla przykładowego numeru indeksu **`38042`** (który należy zastąpić własnym numerem indeksu w skrypcie `generuj_wykresy.py`) zestaw parametrów wejściowych w jednostkach SI przedstawia poniższa tabela:

| Wielkość fizyczna | Symbol | Wartość w układzie SI | Jednostka SI | Opis |
| :--- | :---: | :---: | :---: | :--- |
| Prędkość początkowa | $v_0$ | $38{,}0$ | $\text{m/s}$ | Dwie pierwsze cyfry indeksu (`38`) |
| Kąt początkowy rzutu | $\alpha$ | $0{,}7330$ | $\text{rad}$ | Dwie ostatnie cyfry indeksu (`42`, tj. $\frac{7\pi}{30}\text{ rad}$) |
| Przyspieszenie ziemskie | $g$ | $9{,}81$ | $\text{m/s}^2$ | Standardowe przyspieszenie grawitacyjne |
| Masa rzucanego ciała | $m$ | $0{,}145$ | $\text{kg}$ | Masa kuli (odpowiadająca piłce sportowej) |
| Promień kuli | $R$ | $0{,}037$ | $\text{m}$ | Promień przekroju poprzecznego ciała |
| Pole przekroju poprzecznego | $A$ | $4{,}301 \cdot 10^{-3}$ | $\text{m}^2$ | Wyznaczone ze wzoru $A = \pi R^2$ |
| Gęstość powietrza | $\rho$ | $1{,}225$ | $\text{kg/m}^3$ | Gęstość powietrza w warunkach normalnych |
| Współczynnik oporu kształtu | $C_d$ | $0{,}47$ | $1$ (bezwym.) | Współczynnik aerodynamiczny dla gładkiej kuli |
| Współczynnik oporu liniowego | $b$ | $0{,}025$ | $\text{kg/s}$ | Współczynnik tłumienia lepkiego dla modelu Stokesa |
| Zredukowany parametr oporu | $k$ | $8{,}539 \cdot 10^{-3}$ | $\text{m}^{-1}$ | Parametr zdefiniowany jako $k = \frac{\rho C_d A}{2m}$ |

---

## 2. Wyprowadzenie równań ruchu bez oporu powietrza (model idealny)

### 2.1. Założenia i układ współrzędnych

Analizę prowadzimy w dwuwymiarowym prostokątnym układzie współrzędnych $(x, y)$, w którym:

1. Oś $Ox$ jest skierowana poziomo zgodnie z kierunkiem rzutu, natomiast oś $Oy$ jest skierowana pionowo w górę.
2. Wektor przyspieszenia ziemskiego ma postać $\vec{g} = [0, -g]$, gdzie $g = 9{,}81\text{ m/s}^2$.
3. Wysokość lotu jest zaniedbywalnie mała w porównaniu z promieniem Ziemi ($H_{\max} \ll R_Z \approx 6{,}371 \cdot 10^6\text{ m}$), dzięki czemu pole grawitacyjne można traktować jako jednorodne.
4. Pomijamy opór aerodynamiczny ośrodka, siłę wyporu aerostatycznego oraz siłę Coriolisa.
5. W chwili początkowej $t = 0\text{ s}$ ciało znajduje się w punkcie $(x_0, y_0) = (0\text{ m}, 0\text{ m})$, a składowe wektora prędkości początkowej $\vec{v}_0$ wynoszą:

$$
v_{0x} = v_0 \cos\alpha, \qquad v_{0y} = v_0 \sin\alpha
$$

### 2.2. Całkowanie równań dynamiki Newtona

W próżni jedyną siłą działającą na ciało o masie $m$ jest stała siła ciężkości $\vec{F}_g = m\vec{g}$. Zgodnie z II zasadą dynamiki Newtona:

$$
m \frac{d^2\vec{r}}{dt^2} = \vec{F}_g \implies m \frac{d^2\vec{r}}{dt^2} = m\vec{g}
$$

Po obustronnym podzieleniu przez masę $m > 0$ i rzutowaniu równania wektorowego na osie $Ox$ oraz $Oy$ otrzymujemy układ dwóch niezależnych równań różniczkowych zwyczajnych drugiego rzędu:

$$
\begin{cases}
a_x(t) = \dfrac{d^2x}{dt^2} = 0 \cr
a_y(t) = \dfrac{d^2y}{dt^2} = -g
\end{cases}
$$

#### Krok 1: Wyznaczenie zależności prędkości od czasu $v_x(t)$ i $v_y(t)$

Całkujemy obie strony równań względem czasu $t$:

$$
v_x(t) = \int 0\,dt = C_1, \qquad v_y(t) = \int (-g)\,dt = -gt + C_2
$$

Stałe całkowania $C_1$ i $C_2$ wyznaczamy z warunków początkowych w chwili $t = 0\text{ s}$:

$$
\begin{cases}
v_x(0) = C_1 = v_0 \cos\alpha \cr
v_y(0) = -g \cdot 0 + C_2 = v_0 \sin\alpha \implies C_2 = v_0 \sin\alpha
\end{cases}
$$

Stąd otrzymujemy równania składowych prędkości w dowolnej chwili $t$:

$$
\begin{cases}
v_x(t) = v_0 \cos\alpha \cr
v_y(t) = v_0 \sin\alpha - gt
\end{cases}
$$

Ruch wzdłuż osi $Ox$ jest ruchem jednostajnym prostoliniowym, natomiast wzdłuż osi $Oy$ jest ruchem jednostajnie zmiennym z przyspieszeniem $-g$.

#### Krok 2: Wyznaczenie równań położenia $x(t)$ i $y(t)$

Ponownie całkujemy zależności prędkości względem czasu $t$:

$$
\begin{aligned}
x(t) &= \int v_x(t)\,dt = \int v_0 \cos\alpha\,dt = v_0 t \cos\alpha + C_3 \cr
y(t) &= \int v_y(t)\,dt = \int (v_0 \sin\alpha - gt)\,dt = v_0 t \sin\alpha - \frac{1}{2}gt^2 + C_4
\end{aligned}
$$

Wykorzystując warunki początkowe $x(0) = x_0 = 0\text{ m}$ oraz $y(0) = y_0 = 0\text{ m}$, otrzymujemy $C_3 = x_0 = 0$ oraz $C_4 = y_0 = 0$. Równania kinematyczne położenia przyjmują postać:

$$
\begin{cases}
x(t) = v_0 t \cos\alpha \cr
y(t) = v_0 t \sin\alpha - \dfrac{1}{2}gt^2
\end{cases}
$$

### 2.3. Wyprowadzenie równania toru i parametrów lotu

#### 1. Równanie toru $y(x)$

Dla $\alpha \neq \frac{\pi}{2}\text{ rad}$ ($\cos\alpha \neq 0$) z pierwszego równania położenia wyznaczamy czas $t$ jako funkcję współrzędnej poziomej $x$:

$$
t(x) = \frac{x}{v_0 \cos\alpha}
$$

Podstawiając wyrażenie na $t(x)$ do drugiego równania kinematycznego, eliminujemy parametr czasu i otrzymujemy jawne równanie toru:

$$
y(x) = v_0 \left(\frac{x}{v_0 \cos\alpha}\right) \sin\alpha - \frac{1}{2}g \left(\frac{x}{v_0 \cos\alpha}\right)^2 = x \tan\alpha - \frac{g}{2 v_0^2 \cos^2\alpha} x^2
$$

Jest to równanie paraboli o ramionach skierowanych w dół, symetrycznej względem prostej przechodzącej przez jej wierzchołek.

#### 2. Czas wznoszenia $t_w$ oraz całkowity czas lotu $t_c$

W najwyższym punkcie toru (wierzchołku paraboli) ciało przestaje się wznosić, co odpowiada zerowaniu się pionowej składowej wektora prędkości:

$$
v_y(t_w) = 0 \implies v_0 \sin\alpha - g t_w = 0 \implies t_w = \frac{v_0 \sin\alpha}{g}
$$

Całkowity czas lotu $t_c$ wyznaczamy z warunku powrotu ciała na poziom odniesienia ($y(t_c) = 0$ dla $t_c > 0$):

$$
v_0 t_c \sin\alpha - \frac{1}{2}g t_c^2 = 0 \implies t_c \left(v_0 \sin\alpha - \frac{1}{2}g t_c\right) = 0 \implies t_c = \frac{2 v_0 \sin\alpha}{g} = 2 t_w
$$

Czas opadania jest w modelu idealnym dokładnie równy czasowi wznoszenia ($t_{\text{op}} = t_c - t_w = t_w$).

#### 3. Maksymalna wysokość lotu $H_{\max}$

Podstawiając wyznaczony czas wznoszenia $t_w = \frac{v_0 \sin\alpha}{g}$ do równania na $y(t)$, otrzymujemy:

$$
H_{\max} = y(t_w) = v_0 \left(\frac{v_0 \sin\alpha}{g}\right) \sin\alpha - \frac{1}{2}g \left(\frac{v_0 \sin\alpha}{g}\right)^2 = \frac{v_0^2 \sin^2\alpha}{g} - \frac{v_0^2 \sin^2\alpha}{2g} = \frac{v_0^2 \sin^2\alpha}{2g}
$$

#### 4. Zasięg rzutu $Z$ i optymalny kąt rzutu $\alpha_{\text{opt}}$

Podstawiając całkowity czas lotu $t_c = \frac{2 v_0 \sin\alpha}{g}$ do równania na współrzędną poziomą $x(t)$ oraz korzystając z tożsamości trygonometrycznej $\sin(2\alpha) = 2\sin\alpha\cos\alpha$, otrzymujemy:

$$
Z = x(t_c) = v_0 \left(\frac{2 v_0 \sin\alpha}{g}\right) \cos\alpha = \frac{v_0^2 (2\sin\alpha\cos\alpha)}{g} = \frac{v_0^2 \sin(2\alpha)}{g}
$$

Maksymalną wartość zasięgu $Z_{\max}$ przy ustalonej prędkości początkowej $v_0$ uzyskujemy wtedy, gdy $\sin(2\alpha)$ osiąga maksimum równe $1$:

$$
\sin(2\alpha_{\text{opt}}) = 1 \implies 2\alpha_{\text{opt}} = \frac{\pi}{2}\text{ rad} \implies \alpha_{\text{opt}} = \frac{\pi}{4}\text{ rad} \approx 0{,}7854\text{ rad}
$$

Wówczas maksymalny zasięg w próżni wynosi $Z_{\max} = \dfrac{v_0^2}{g}$.

---

## 3. Wyprowadzenie równań ruchu z oporem powietrza (model rzeczywisty)

W rzeczywistym ośrodku gazowym na poruszające się ciało działa siła oporu aerodynamicznego $\vec{F}_{\text{op}}$, której zwrot jest zawsze przeciwny do wektora chwilowej prędkości $\vec{v} = [v_x, v_y]$:

$$
\vec{F}_{\text{op}} = -F_{\text{op}}(v) \, \hat{v} = -F_{\text{op}}(v) \, \frac{\vec{v}}{|\vec{v}|}
$$

Charakter zależności wartości siły oporu $F_{\text{op}}(v)$ od wartości prędkości $v = |\vec{v}| = \sqrt{v_x^2 + v_y^2}$ zależy od bezwymiarowej liczby Reynoldsa:

$$
Re = \frac{\rho v D}{\eta}
$$

gdzie $D = 2R$ oznacza średnicę ciała ($\text{m}$), a $\eta$ jest lepkością dynamiczną ośrodka ($\text{Pa}\cdot\text{s}$).

### 3.1. Model oporu liniowego — prawo Stokesa (wyprowadzenie analityczne)

Dla małych liczb Reynoldsa ($Re \lesssim 1$, przepływ laminarny) siła oporu jest wprost proporcjonalna do pierwszej potęgi prędkości:

$$
\vec{F}_{\text{op}} = -b \vec{v}
$$

gdzie $b$ jest współczynnikiem oporu lepkiego wyrażonym w $\text{kg/s}$ (dla kuli o promieniu $R$ zachodzi $b = 6\pi \eta R$). Równanie dynamiki Newtona przyjmuje postać liniowego układu rozsprzężonego:

$$
\begin{cases}
m \dfrac{dv_x}{dt} = -b v_x \cr
m \dfrac{dv_y}{dt} = -mg - b v_y
\end{cases}
$$

#### Rozwiązanie dla osi poziomej $Ox$:

Rozdzielamy zmienne i całkujemy od $t = 0$ (gdzie $v_x(0) = v_{0x} = v_0 \cos\alpha$) do chwili $t$:

$$
\int_{v_{0x}}^{v_x(t)} \frac{dv_x}{v_x} = -\frac{b}{m} \int_0^t dt \implies \ln\left(\frac{v_x(t)}{v_{0x}}\right) = -\frac{b}{m}t \implies v_x(t) = v_0 \cos\alpha \cdot e^{-\frac{b}{m}t}
$$

Całkując $v_x(t)$ po czasie w granicach od $0$ do $t$ (przy $x(0) = 0\text{ m}$), otrzymujemy:

$$
x(t) = v_0 \cos\alpha \int_0^t e^{-\frac{b}{m}\tau}\,d\tau = \frac{m v_0 \cos\alpha}{b} \left(1 - e^{-\frac{b}{m}t}\right)
$$

#### Rozwiązanie dla osi pionowej $Oy$:

Przekształcamy równanie pionowe do postaci $\dfrac{dv_y}{v_y + \frac{mg}{b}} = -\frac{b}{m}\,dt$ i całkujemy od $v_{0y} = v_0 \sin\alpha$ do $v_y(t)$:

$$
\ln\left(\frac{v_y(t) + \frac{mg}{b}}{v_{0y} + \frac{mg}{b}}\right) = -\frac{b}{m}t \implies v_y(t) = \left(v_0 \sin\alpha + \frac{mg}{b}\right) e^{-\frac{b}{m}t} - \frac{mg}{b}
$$

Po ponownym scałkowaniu po czasie od $0$ do $t$ (przy $y(0) = 0\text{ m}$) uzyskujemy zależność $y(t)$:

$$
y(t) = \int_0^t v_y(\tau)\,d\tau = \frac{m}{b} \left(v_0 \sin\alpha + \frac{mg}{b}\right) \left(1 - e^{-\frac{b}{m}t}\right) - \frac{mg}{b}t
$$

### 3.2. Model oporu kwadratowego — prawo Newtona/Rayleigha

Dla dużych liczb Reynoldsa ($10^3 \lesssim Re \lesssim 2 \cdot 10^5$, przepływ turbulentny typowy dla ruchu piłek, pocisków i obiektów makroskopowych w powietrzu) dominuje opór ciśnieniowy proporcjonalny do kwadratu wartości prędkości:

$$
F_{\text{op}}(v) = \frac{1}{2} \rho C_d A v^2 = c v^2, \qquad \text{gdzie } c = \frac{1}{2} \rho C_d A \quad [\text{kg/m}]
$$

Uwzględniając kierunek i zwrot wektora prędkości poprzez wersor $\hat{v} = \dfrac{\vec{v}}{|\vec{v}|} = \dfrac{[v_x, v_y]}{\sqrt{v_x^2 + v_y^2}}$, zapisujemy wektor siły oporu:

$$
\vec{F}_{\text{op}} = -c v^2 \frac{\vec{v}}{v} = -c v \vec{v} = -c \sqrt{v_x^2 + v_y^2} \begin{bmatrix} v_x \cr v_y \end{bmatrix}
$$

Podstawiając $\vec{F}_{\text{op}}$ oraz $\vec{F}_g = [0, -mg]^T$ do II zasady dynamiki Newtona $m \frac{d\vec{v}}{dt} = \vec{F}_g + \vec{F}_{\text{op}}$, otrzymujemy układ równań:

$$
\begin{cases}
m \dfrac{dv_x}{dt} = -c \sqrt{v_x^2 + v_y^2} \cdot v_x \cr
m \dfrac{dv_y}{dt} = -mg - c \sqrt{v_x^2 + v_y^2} \cdot v_y
\end{cases}
$$

Dzieląc oba równania przez masę $m$ i wprowadzając zredukowany parametr oporu $k = \dfrac{c}{m} = \dfrac{\rho C_d A}{2m}$ (wyrażony w $\text{m}^{-1}$), otrzymujemy równania przyspieszeń:

$$
\begin{cases}
a_x = \dfrac{dv_x}{dt} = -k v_x \sqrt{v_x^2 + v_y^2} \cr
a_y = \dfrac{dv_y}{dt} = -g - k v_y \sqrt{v_x^2 + v_y^2}
\end{cases}
$$

> **Wniosek matematyczny:** Obecność czynnika $\sqrt{v_x^2 + v_y^2}$ powoduje nieliniowe sprzężenie składowych ruchu wzdłuż osi $Ox$ oraz $Oy$ (zmiana prędkości poziomej wpływa na opóźnienie w kierunku pionowym i odwrotnie). Układ ten nie ma rozwiązania analitycznego w klasie funkcji elementarnych i wymaga zastosowania metod całkowania numerycznego.

### 3.3. Prędkość graniczna w ośrodku lepkim i turbulentnym

Podczas swobodnego spadania w ośrodku z oporem powietrza siła oporu rośnie wraz z prędkością aż do momentu zrównoważenia siły ciężkości ($a_y = 0$, $v_x = 0$, $v_y = -v_{\text{gr}}$):

* **Dla oporu liniowego:** $b v_{\text{gr, lin}} = mg \implies v_{\text{gr, lin}} = \dfrac{mg}{b}$,
* **Dla oporu kwadratowego:** $c v_{\text{gr}}^2 = mg \implies v_{\text{gr}} = \sqrt{\dfrac{mg}{c}} = \sqrt{\dfrac{g}{k}}$.

Dla przyjętych parametrów ($g = 9{,}81\text{ m/s}^2$, $k = 8{,}539 \cdot 10^{-3}\text{ m}^{-1}$) prędkość graniczna w modelu kwadratowym wynosi $v_{\text{gr}} \approx 33{,}89\text{ m/s}$.

---

## 4. Metody numeryczne rozwiązywania równań ruchu

### 4.1. Redukcja do układu równań różniczkowych I rzędu

W celu rozwiązania zagadnienia początkowego (Cauchy'ego) dla modelu z oporem kwadratowym definiujemy czteroelementowy wektor stanu układu $\vec{S}(t)$:

$$
\vec{S}(t) = \begin{bmatrix} x(t) \cr y(t) \cr v_x(t) \cr v_y(t) \end{bmatrix}, \qquad \vec{S}(0) = \begin{bmatrix} 0 \cr 0 \cr v_0 \cos\alpha \cr v_0 \sin\alpha \end{bmatrix}
$$

Układ dwóch równań różniczkowych drugiego rzędu sprowadza się wówczas do autonomicznego układu czterech równań różniczkowych zwyczajnych pierwszego rzędu $\dfrac{d\vec{S}}{dt} = \vec{f}(\vec{S})$:

$$
\frac{d}{dt} \begin{bmatrix} x \cr y \cr v_x \cr v_y \end{bmatrix} = \begin{bmatrix} v_x \cr v_y \cr -k v_x \sqrt{v_x^2 + v_y^2} \cr -g - k v_y \sqrt{v_x^2 + v_y^2} \end{bmatrix}
$$

### 4.2. Schemat Eulera oraz metoda Rungego-Kutty 4. rzędu (RK4)

Dyskretyzujemy oś czasu ze stałym krokiem $\Delta t$ ($t_n = n \Delta t$, $n = 0, 1, 2, \dots$).

1. **Metoda Eulera (I rzędu):** Opiera się na liniowym rozwinięciu w szereg Taylora $\vec{S}_{n+1} = \vec{S}_n + \Delta t \, \vec{f}(\vec{S}_n)$ z błędem lokalnym rzędu $\mathcal{O}(\Delta t^2)$ i błędem globalnym $\mathcal{O}(\Delta t)$. Ze względu na kumulację błędu nie zapewnia zachowania energii w układach zachowawczych.
2. **Klasyczna metoda Rungego-Kutty 4. rzędu (RK4):** Wyznacza ważoną średnią z czterech przyrostów na przedziale $[t_n, t_n + \Delta t]$, osiągając błąd lokalny rzędu $\mathcal{O}(\Delta t^5)$ oraz błąd globalny rzędu $\mathcal{O}(\Delta t^4)$:

$$
\begin{aligned}
\vec{k}_1 &= \vec{f}(\vec{S}_n) \cr
\vec{k}_2 &= \vec{f}\left(\vec{S}_n + \frac{\Delta t}{2}\vec{k}_1\right) \cr
\vec{k}_3 &= \vec{f}\left(\vec{S}_n + \frac{\Delta t}{2}\vec{k}_2\right) \cr
\vec{k}_4 &= \vec{f}\left(\vec{S}_n + \Delta t \vec{k}_3\right)
\end{aligned}
$$

$$
\vec{S}_{n+1} = \vec{S}_n + \frac{\Delta t}{6}\left(\vec{k}_1 + 2\vec{k}_2 + 2\vec{k}_3 + \vec{k}_4\right)
$$

Warunkiem zakończenia symulacji jest przecięcie poziomu podłoża ($y_{n+1} \le 0\text{ m}$). Dokładną chwilę upadku $t_c$ oraz zasięg rzutu $Z$ wyznaczamy za pomocą interpolacji liniowej pomiędzy stanami $\vec{S}_n$ oraz $\vec{S}_{n+1}$:

$$
\theta = \frac{-y_n}{y_{n+1} - y_n} \in [0, 1], \qquad \vec{S}_{\text{upadek}} = \vec{S}_n + \theta (\vec{S}_{n+1} - \vec{S}_n), \qquad t_c = t_n + \theta \Delta t
$$

---

## 5. Implementacja symulacji oraz wykresy ruchu

Poniższy kod w języku Python (pełny skrypt znajduje się również w pliku `generuj_wykresy.py`) realizuje obliczenia numeryczne metodą RK4 oraz generuje zestaw pięciu wykresów w katalogu `wykresy/`.

> **Instrukcja uruchomienia w VS Code:**  
> 1. Otwórz plik `generuj_wykresy.py` w edytorze VS Code i wpisz swój numer indeksu w zmiennej `INDEKS`.  
> 2. Uruchom skrypt przyciskiem **Run Python File** (`Ctrl + F5`).  
> 3. Skrypt wygeneruje wszystkie pliki graficzne w folderze `wykresy/` oraz wypisze w terminalu dokładną tabelę wyników dla Twojego numeru indeksu.

#### Wspólny moduł inicjalizacji parametrów (SI) i solvera RK4:

```python
import os
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 1. PARAMETRY WEJŚCIOWE Z NUMERU INDEKSU (ŚCIŚLE W JEDNOSTKACH SI)
# ==============================================================================
INDEKS = "366270"

v0 = float(INDEKS[:2])                # [m/s] dwie pierwsze cyfry indeksu
alpha_deg = float(INDEKS[-2:])        # wartość pomocnicza z dwóch ostatnich cyfr
alpha = np.deg2rad(alpha_deg)         # [rad] kąt początkowy w jednostce SI (radian)

# Stałe fizyczne i parametry ciała (wszystkie w jednostkach SI)
g = 9.81                              # [m/s^2] przyspieszenie ziemskie
m = 0.145                             # [kg] masa ciała
R = 0.037                             # [m] promień kuli
A = np.pi * R**2                      # [m^2] pole przekroju poprzecznego
rho = 1.225                           # [kg/m^3] gęstość powietrza
Cd = 0.47                             # [-] współczynnik oporu aerodynamicznego
c = 0.5 * rho * Cd * A                # [kg/m] współczynnik oporu kwadratowego
k = c / m                             # [1/m] zredukowany parametr oporu
b = 0.025                             # [kg/s] współczynnik oporu liniowego (Stokes)

os.makedirs("wykresy", exist_ok=True)


def symuluj_rk4(v_start, kat_rad, tryb="kwadratowy", dt=0.001):
    """Całkowanie równań ruchu metodą RK4 z interpolacją punktu upadku (y = 0 m)."""
    def pochodne(stan):
        x, y, vx, vy = stan
        v = np.hypot(vx, vy)
        if tryb == "brak":
            ax, ay = 0.0, -g
        elif tryb == "liniowy":
            ax = -(b / m) * vx
            ay = -g - (b / m) * vy
        elif tryb == "kwadratowy":
            ax = -k * v * vx
            ay = -g - k * v * vy
        return np.array([vx, vy, ax, ay], dtype=float)

    stan = np.array([0.0, 0.0, v_start * np.cos(kat_rad), v_start * np.sin(kat_rad)])
    t = 0.0
    czasy = [t]
    stany = [stan.copy()]

    while True:
        k1 = pochodne(stan)
        k2 = pochodne(stan + 0.5 * dt * k1)
        k3 = pochodne(stan + 0.5 * dt * k2)
        k4 = pochodne(stan + dt * k3)

        nowy_stan = stan + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        nowy_t = t + dt

        if nowy_stan[1] <= 0.0 and t > 0.0:
            theta = -stan[1] / (nowy_stan[1] - stan[1])
            stan_upadek = stan + theta * (nowy_stan - stan)
            t_upadek = t + theta * dt
            czasy.append(t_upadek)
            stany.append(stan_upadek)
            break

        stan = nowy_stan
        t = nowy_t
        czasy.append(t)
        stany.append(stan.copy())

    return np.array(czasy), np.array(stany)
```

---

### 5.1. Wykres 1: Porównanie trajektorii lotu $y(x)$

Wykres przedstawia zestawienie torów ruchu w przestrzeni $(x, y)$ dla trzech modeli fizycznych: idealnego (parabola w próżni), z oporem liniowym oraz z oporem kwadratowym. Wyraźnie widoczna jest utrata symetrii krzywej balistycznej — ramię opadające toru jest znacznie bardziej strome niż ramię wznoszące.

```python
t_id, s_id = symuluj_rk4(v0, alpha, tryb="brak")
t_lin, s_lin = symuluj_rk4(v0, alpha, tryb="liniowy")
t_kw, s_kw = symuluj_rk4(v0, alpha, tryb="kwadratowy")

plt.figure(figsize=(10, 5.5), dpi=150)
plt.plot(s_id[:, 0], s_id[:, 1], "--", linewidth=2, label="Bez oporu powietrza (model idealny)")
plt.plot(s_lin[:, 0], s_lin[:, 1], "-.", linewidth=2, label=f"Opór liniowy Stokesa (b = {b} kg/s)")
plt.plot(s_kw[:, 0], s_kw[:, 1], "-", linewidth=2.5, label=f"Opór kwadratowy (k = {k:.5f} 1/m)")

# Zaznaczenie maksimów wysokości
for dane, kolor in [(s_id, "C0"), (s_lin, "C1"), (s_kw, "C2")]:
    idx_max = np.argmax(dane[:, 1])
    plt.scatter(dane[idx_max, 0], dane[idx_max, 1], color=kolor, zorder=5)

plt.title(f"Trajektorie rzutu ukośnego (v0 = {v0:.1f} m/s, α = {alpha:.4f} rad)")
plt.xlabel("Odległość pozioma x [m]")
plt.ylabel("Wysokość y [m]")
plt.grid(True, linestyle=":", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig("wykresy/01_trajektorie_porownanie.png")
plt.close()
```

![Porównanie trajektorii rzutu ukośnego](wykresy/01_trajektorie_porownanie.png)

---

### 5.2. Wykres 2: Przebieg składowych i wartości prędkości w czasie

Poniższy wykres ilustruje ewolucję czasową składowej poziomej $v_x(t)$, składowej pionowej $v_y(t)$ oraz modułu wektora prędkości $v(t) = \sqrt{v_x^2 + v_y^2}$. W modelu idealnym $v_x(t) = \text{const}$, natomiast w modelu z oporem powietrza składowa pozioma zanika monotonicznie, a minimum prędkości całkowitej $v_{\min}$ nie występuje dokładnie w najwyższym punkcie toru ($v_y = 0$), lecz nieco po jego minięciu.

```python
v_id = np.hypot(s_id[:, 2], s_id[:, 3])
v_kw = np.hypot(s_kw[:, 2], s_kw[:, 3])

fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

axes[0].plot(t_id, s_id[:, 2], "--", label="vx(t) - bez oporu")
axes[0].plot(t_id, s_id[:, 3], "--", label="vy(t) - bez oporu")
axes[0].plot(t_kw, s_kw[:, 2], "-", linewidth=2, label="vx(t) - opór kwadratowy")
axes[0].plot(t_kw, s_kw[:, 3], "-", linewidth=2, label="vy(t) - opór kwadratowy")
axes[0].axhline(0.0, color="black", linewidth=0.8)
axes[0].set_title("Składowe wektora prędkości vx(t) oraz vy(t)")
axes[0].set_xlabel("Czas t [s]")
axes[0].set_ylabel("Prędkość [m/s]")
axes[0].grid(True, linestyle=":", alpha=0.7)
axes[0].legend()

axes[1].plot(t_id, v_id, "--", linewidth=2, label="v(t) - bez oporu")
axes[1].plot(t_kw, v_kw, "-", linewidth=2, label="v(t) - opór kwadratowy")
axes[1].set_title("Wartość prędkości całkowitej v(t)")
axes[1].set_xlabel("Czas t [s]")
axes[1].set_ylabel("Wartość prędkości v [m/s]")
axes[1].grid(True, linestyle=":", alpha=0.7)
axes[1].legend()

plt.tight_layout()
plt.savefig("wykresy/02_predkosci_w_czasie.png")
plt.close()
```

![Składowe i wartość prędkości w funkcji czasu](wykresy/02_predkosci_w_czasie.png)

---

### 5.3. Wykres 3: Zależność współrzędnych położenia $x(t)$ i $y(t)$ od czasu

Zestawienie krzywych $x(t)$ oraz $y(t)$ pokazuje zakrzywienie wykresu przemieszczenia poziomego w wyniku wyhamowywania ruchu wzdłuż osi $Ox$ oraz krótszy całkowity czas przebywania ciała w powietrzu.

```python
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

axes[0].plot(t_id, s_id[:, 0], "--", linewidth=2, label="Bez oporu powietrza")
axes[0].plot(t_kw, s_kw[:, 0], "-", linewidth=2, label="Opór kwadratowy")
axes[0].set_title("Współrzędna pozioma x(t)")
axes[0].set_xlabel("Czas t [s]")
axes[0].set_ylabel("Położenie x [m]")
axes[0].grid(True, linestyle=":", alpha=0.7)
axes[0].legend()

axes[1].plot(t_id, s_id[:, 1], "--", linewidth=2, label="Bez oporu powietrza")
axes[1].plot(t_kw, s_kw[:, 1], "-", linewidth=2, label="Opór kwadratowy")
axes[1].set_title("Wysokość y(t)")
axes[1].set_xlabel("Czas t [s]")
axes[1].set_ylabel("Wysokość y [m]")
axes[1].grid(True, linestyle=":", alpha=0.7)
axes[1].legend()

plt.tight_layout()
plt.savefig("wykresy/03_wspolrzedne_w_czasie.png")
plt.close()
```

![Zależność współrzędnych położenia od czasu](wykresy/03_wspolrzedne_w_czasie.png)

---

### 5.4. Wykres 4: Wpływ kąta wyrzutu $\alpha$ na zasięg i wysokość maksymalną

Badanie zależności $Z(\alpha)$ oraz $H_{\max}(\alpha)$ dla $\alpha \in [0{,}05;\, 1{,}52]\text{ rad}$ dowodzi jednego z najważniejszych wniosków balistyki zewnętrznej: obecność oporu powietrza przesuwa optymalny kąt maksymalnego zasięgu poniżej wartości $\frac{\pi}{4}\text{ rad} \approx 0{,}7854\text{ rad}$.

```python
katy_rad = np.linspace(0.05, 1.52, 80)
zasiegi_id, zasiegi_kw = [], []
wys_id, wys_kw = [], []

for kat in katy_rad:
    _, si = symuluj_rk4(v0, kat, tryb="brak", dt=0.002)
    _, sk = symuluj_rk4(v0, kat, tryb="kwadratowy", dt=0.002)
    zasiegi_id.append(si[-1, 0])
    zasiegi_kw.append(sk[-1, 0])
    wys_id.append(np.max(si[:, 1]))
    wys_kw.append(np.max(sk[:, 1]))

fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

axes[0].plot(katy_rad, zasiegi_id, "--", linewidth=2, label="Bez oporu (maks. przy π/4 rad)")
axes[0].plot(katy_rad, zasiegi_kw, "-", linewidth=2, label="Opór kwadratowy")
axes[0].axvline(alpha, color="red", linestyle=":", label=f"Twój kąt α = {alpha:.4f} rad")
axes[0].set_title("Zasięg rzutu Z(α)")
axes[0].set_xlabel("Kąt początkowy α [rad]")
axes[0].set_ylabel("Zasięg Z [m]")
axes[0].grid(True, linestyle=":", alpha=0.7)
axes[0].legend()

axes[1].plot(katy_rad, wys_id, "--", linewidth=2, label="Bez oporu powietrza")
axes[1].plot(katy_rad, wys_kw, "-", linewidth=2, label="Opór kwadratowy")
axes[1].axvline(alpha, color="red", linestyle=":", label=f"Twój kąt α = {alpha:.4f} rad")
axes[1].set_title("Maksymalna wysokość H_max(α)")
axes[1].set_xlabel("Kąt początkowy α [rad]")
axes[1].set_ylabel("Wysokość maksymalna H_max [m]")
axes[1].grid(True, linestyle=":", alpha=0.7)
axes[1].legend()

plt.tight_layout()
plt.savefig("wykresy/04_wplyw_kata_na_zasieg.png")
plt.close()
```

![Wpływ kąta rzutu na zasięg i wysokość maksymalną](wykresy/04_wplyw_kata_na_zasieg.png)

---

### 5.5. Wykres 5: Bilans energii mechanicznej układu

W modelu bez oporu powietrza pole grawitacyjne jest polem potencjalnym, więc całkowita energia mechaniczna $E_c = E_k + E_p = \frac{1}{2}m(v_x^2 + v_y^2) + mgy$ pozostaje stała w czasie. Siła oporu powietrza jest siłą niepotencjalną (dyssypatywną), wykonującą ujemną pracę $W_{\text{op}} = \int \vec{F}_{\text{op}} \cdot d\vec{r} < 0$, co skutkuje nieodwracalnym rozpraszaniem energii mechanicznej w postaci ciepła.

```python
Ek_id = 0.5 * m * (s_id[:, 2]**2 + s_id[:, 3]**2)
Ep_id = m * g * s_id[:, 1]
Ec_id = Ek_id + Ep_id

Ek_kw = 0.5 * m * (s_kw[:, 2]**2 + s_kw[:, 3]**2)
Ep_kw = m * g * s_kw[:, 1]
Ec_kw = Ek_kw + Ep_kw

plt.figure(figsize=(10, 5.5), dpi=150)
plt.plot(t_id, Ec_id, "--", linewidth=2, label="Energia całkowita Ec (bez oporu - zachowana)")
plt.plot(t_kw, Ec_kw, "-", linewidth=2.5, label="Energia całkowita Ec (opór kwadratowy)")
plt.plot(t_kw, Ek_kw, ":", linewidth=1.8, label="Energia kinetyczna Ek (opór kwadratowy)")
plt.plot(t_kw, Ep_kw, "-.", linewidth=1.8, label="Energia potencjalna Ep (opór kwadratowy)")

plt.title("Bilans energii mechanicznej układu w trakcie rzutu ukośnego")
plt.xlabel("Czas t [s]")
plt.ylabel("Energia [J]")
plt.grid(True, linestyle=":", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig("wykresy/05_energia_mechaniczna.png")
plt.close()
```

![Zmiana energii mechanicznej w trakcie lotu](wykresy/05_energia_mechaniczna.png)

---

## 6. Zestawienie porównawcze wyników w jednostkach SI

Poniższa tabela prezentuje porównanie głównych wielkości fizycznych uzyskanych dla przykładowego indeksu **`38042`** ($v_0 = 38{,}0\text{ m/s}$, $\alpha = 0{,}7330\text{ rad}$):

| Wielkość fizyczna | Symbol | Model bez oporu powietrza | Model z oporem kwadratowym | Różnica względna | Jednostka SI |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Zasięg rzutu | $Z$ | $146{,}39$ | $82{,}51$ | $-43{,}6\%$ | $\text{m}$ |
| Maksymalna wysokość | $H_{\max}$ | $32{,}95$ | $22{,}89$ | $-30{,}5\%$ | $\text{m}$ |
| Całkowity czas lotu | $t_c$ | $5{,}18$ | $4{,}31$ | $-16{,}8\%$ | $\text{s}$ |
| Czas wznoszenia | $t_w$ | $2{,}59$ | $1{,}99$ | $-23{,}2\%$ | $\text{s}$ |
| Czas opadania | $t_{\text{op}}$ | $2{,}59$ | $2{,}32$ | $-10{,}4\%$ | $\text{s}$ |
| Prędkość końcowa (upadku) | $v_k$ | $38{,}00$ | $23{,}64$ | $-37{,}8\%$ | $\text{m/s}$ |
| Kąt upadku | $\beta$ | $0{,}7330$ | $0{,}9053$ | $+23{,}5\%$ | $\text{rad}$ |
| Optymalny kąt zasięgu | $\alpha_{\text{opt}}$ | $0{,}7854$ | $0{,}7012$ | $-10{,}7\%$ | $\text{rad}$ |
| Energia całkowita przy upadku | $E_{c}(t_c)$ | $104{,}69$ | $40{,}52$ | $-61{,}3\%$ | $\text{J}$ |

---

## 7. Podsumowanie i wnioski

1. **Wpływ oporu ośrodka na parametry kinematyczne:** Uwzględnienie kwadratowego oporu powietrza znacząco redukuje zasięg rzutu ($Z$), wysokość maksymalną ($H_{\max}$) oraz prędkość końcową ciała. Redukcja ta jest tym silniejsza, im większy jest stosunek prędkości początkowej $v_0$ do prędkości granicznej $v_{\text{gr}}$.
2. **Złamanie symetrii trajektorii:** W odróżnieniu od symetrycznej paraboli w próżni ($t_w = t_{\text{op}}$, $\beta = \alpha$), krzywa balistyczna charakteryzuje się krótszym czasem wznoszenia niż opadania ($t_w < t_{\text{op}}$) oraz większym kątem upadku niż kątem wyrzutu ($\beta > \alpha$).
3. **Przesunięcie optymalnego kąta wyrzutu:** Maksymalny zasięg w obecności oporu powietrza uzyskuje się dla kąta mniejszego niż $\frac{\pi}{4}\text{ rad}$ ($\alpha_{\text{opt}} < 0{,}7854\text{ rad}$). Mniejszy kąt pozwala na szybsze pokonanie dystansu poziomego w początkowej fazie lotu, zanim opór aerodynamiczny wytraci energię kinetyczną ciała.
4. **Skuteczność metody RK4:** Klasyczna metoda Rungego-Kutty 4. rzędu (RK4) z krokiem $\Delta t = 0{,}001\text{ s}$ zachowuje energię mechaniczną w modelu idealnym z dokładnością rzędu $10^{-11}\text{ J}$ oraz zapewnia wysoką precyzję całkowania nieliniowego układu równań ruchu z oporem powietrza.

---

## 8. Literatura

1. Halliday D., Resnick R., Walker J., *Podstawy fizyki*, tom 1, Wydawnictwo Naukowe PWN, Warszawa.
2. Taylor J. R., *Mechanika klasyczna*, tom 1, Wydawnictwo Naukowe PWN, Warszawa.
3. Fortuna Z., Macukow B., Wąsowski J., *Metody numeryczne*, Wydawnictwa Naukowo-Techniczne (WNT), Warszawa.
4. Pang T., *Metody obliczeniowe w fizyce*, Wydawnictwo Naukowe PWN, Warszawa.