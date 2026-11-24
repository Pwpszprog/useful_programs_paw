# Rzut ukośny — od modelu idealnego do ruchu z oporem powietrza



**Temat:** dwuwymiarowy ruch ciała wyrzuconego ukośnie w jednorodnym polu grawitacyjnym bez oporów ruchu oraz z uwzględnieniem oporu powietrza.  
**Metoda:** rozwiązania analityczne, stworzenie symulacji w Pythonie i interpretacja wykresów.  
**Charakter danych:** przykład modelowy; wszystkie wykresy przedstawiają obliczenia symulacyjne, a nie pomiary laboratoryjne.

Galileusz stwierdził, że
> ruch rzuconego ciała można rozłożyć na dwa niezależne ruchy: jednostajny w poziomie i jednostajnie przyspieszony w pionie.

Newton stwierdził, że
> zmiana pędu ciała jest proporcjonalna do przyłożonej siły wypadkowej i odbywa się w kierunku jej działania.

Połączenie tych zasad pozwala opisać tor lotu w próżni za pomocą klasycznej paraboli. W rzeczywistych warunkach pojawia się jednak siła oporu ośrodka, która zaburza symetrię toru, skraca zasięg i wymaga zastosowania metod numerycznych. W tym raporcie porównujemy oba podejścia, analizując położenie, prędkość oraz energię układu.

Wektor położenia mówi, gdzie znajduje się ciało w danej chwili. Wektor prędkości określa, dokąd i jak szybko zmierza. 
Dopiero razem pozwalają jasno określić tor ruchu i zachowanie rzuconego ciała.


## Spis treści

1. [Cel i model fizyczny](#1-cel-i-model-fizyczny)
2. [Rzut ukośny bez oporu powietrza — rozwiązanie analityczne](#2-rzut-ukośny-bez-oporu-powietrza--rozwiązanie-analityczne)
3. [Przykład liczbowy i przebiegi w modelu idealnym](#3-przykład-liczbowy-i-przebiegi-w-modelu-idealnym)
4. [Energia mechaniczna w rzucie idealnym](#4-energia-mechaniczna-w-rzucie-idealnym)
5. [Rzut ukośny z oporem powietrza](#5-rzut-ukośny-z-oporem-powietrza)
6. [Porównanie obu modeli i analiza zasięgu](#6-porównanie-obu-modeli-i-analiza-zasięgu)
7. [Kod i odtwarzanie wykresów](#7-kod-i-odtwarzanie-wykresów)
8. [Jak wykonać doświadczenie](#8-jak-wykonać-doświadczenie)
9. [Wnioski](#9-wnioski)
10. [Dalsza lektura](#10-dalsza-lektura)

---

## 1. Cel i model fizyczny

Dlaczego kamień rzucony pod kątem w próżni zakreśla symetryczną parabolę, a piłka lub pocisk w powietrzu opada znacznie szybciej i bardziej stromo, niż wznosił się na początku? W idealnym modelu jedyną siłą działającą na ciało jest stała siła grawitacji. W rzeczywistości ośrodek stawia opór skierowany zawsze przeciwnie do aktualnego wektora prędkości, co sprzęga ze sobą ruch poziomy i pionowy.

Celem raportu jest wyprowadzenie równań ruchu dla rzutu ukośnego, zbadanie zmian położenia, prędkości i energii w czasie oraz ilościowe porównanie trajektorii idealnej z trajektorią w obecności oporu powietrza.

### 1.1. Założenia

- Ruch odbywa się w pionowej płaszczyźnie $(x, y)$, gdzie oś $x$ skierowana jest poziomo, a oś $y$ pionowo w górę.
- Ciało traktujemy jako punkt materialny o stałej masie $m$.
- Przyspieszenie ziemskie $\vec{g}$ jest stałe co do wartości i kierunku w całym obszarze ruchu (pomijamy krzywiznę Ziemi oraz zmiany wysokości).
- Wyrzut następuje z początku układu współrzędnych $(x_0 = 0,\; y_0 = 0)$ z prędkością początkową $v_0$ pod kątem $\alpha$ do poziomu.
- W modelu idealnym pomijamy oddziaływanie z powietrzem (brak oporu i siły wyporu).
- W modelu z oporem uwzględniamy siłę oporu aerodynamicznego skierowaną przeciwnie do wektora prędkości $\vec{v}$ (rozważamy opór liniowy $\vec{F}_{\mathrm{op}} = -b\vec{v}$ oraz kwadratowy $\vec{F}_{\mathrm{op}} = -c\,v\,\vec{v}$).

Współrzędna $y$ oznacza **wysokość nad poziomem wyrzutu**. Lot kończy się w chwili ponownego osiągnięcia poziomu $y = 0$.

![Schemat rzutu ukośnego: układ współrzędnych, wektor prędkości początkowej pod kątem alfa, tor lotu oraz wektory siły ciężkości i siły oporu]( )

*Rysunek 1. Schemat układu fizycznego dla rzutu ukośnego. W modelu idealnym działa wyłącznie siła ciężkości $\vec{F}_g$; w modelu rzeczywistym pojawia się dodatkowo siła oporu $\vec{F}_{\mathrm{op}}$ przeciwna do wektora prędkości $\vec{v}$. Schemat nie zachowuje skali.*

### 1.2. Symbole i jednostki

Wszystkie wielkości fizyczne w raporcie wyrażono w podstawowych lub pochodnych jednostkach układu **SI**.

| Symbol | Znaczenie | Jednostka SI |
| :--- | :--- | :---: |
| $m$ | masa rzucanego ciała | $\mathrm{kg}$ |
| $g$ | przyspieszenie ziemskie ($g = 9{,}81\,\mathrm{m/s^2}$) | $\mathrm{m/s^2}$ |
| $x$, $y$ | współrzędna pozioma (odległość) i pionowa (wysokość) | $\mathrm{m}$ |
| $v_0$ | wartość prędkości początkowej | $\mathrm{m/s}$ |
| $\alpha$ | kąt rzutu względem poziomu | $\mathrm{rad}$ (oraz $\mathrm{^\circ}$) |
| $v_x$, $v_y$, $v$ | składowe prędkości oraz jej wartość całkowita | $\mathrm{m/s}$ |
| $a_x$, $a_y$ | składowe przyspieszenia | $\mathrm{m/s^2}$ |
| $t$, $t_c$ | czas ruchu oraz całkowity czas lotu | $\mathrm{s}$ |
| $Z$, $H_{\max}$ | zasięg rzutu oraz maksymalna wysokość lotu | $\mathrm{m}$ |
| $E_k$, $E_p$, $E$ | energia kinetyczna, potencjalna i mechaniczna całkowita | $\mathrm{J}$ |
| $b$ | współczynnik oporu liniowego (lepkiego) | $\mathrm{N\,s/m}$ ($\mathrm{kg/s}$) |
| $c$ | współczynnik oporu kwadratowego ($c = \frac{1}{2}\rho C_d S$) | $\mathrm{N\,s^2/m^2}$ ($\mathrm{kg/m}$) |

Nie należy mylić stopni z radianami w obliczeniach numerycznych: funkcje trygonometryczne w równaniach oraz bibliotekach programistycznych wymagają podania kąta $\alpha$ w radianach ($\alpha_{\mathrm{rad}} = \alpha_{\mathrm{deg}} \cdot \frac{\pi}{180^\circ}$).

- Ruch odbywa się w pionowej płaszczyźnie $(x, y)$, gdzie oś $x$ skierowana jest poziomo, a oś $y$ pionowo w górę.
- Ciało traktujemy jako punkt materialny o stałej masie $m$.
-