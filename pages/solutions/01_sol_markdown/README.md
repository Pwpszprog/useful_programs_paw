# Tytuł projektu

## Pierwsza podsekcja

To jest krótki tekst umieszczony w pierwszej podsekcji. Zawiera on **pogrubiony fragment**, *pochylony fragment* oraz ~~przekreślony fragment~~

Przydatne wzory matematyczne w tekście:
- Wzór na energię: $E = mc^2$
- Twierdzenie Pitagorasa: $a^2 + b^2 = c^2$
- Funkcja kwadratowa: $f(x) = x^2$
## Druga podsekcja

To jest kolejny krótki tekst w ramach drugiej podsekcji. Znajdziesz tutaj link do platformy [Google Colab](http://colab.research.google.com)

### Lista punktowana
* Pierwszy element listy
* Drugi element listy
* Trzeci element listy

### Lista numerowana
1. Pierwszy krok
2. Drugi krok
3. Trzeci krok

### Checklista
- [x] Zadanie wykonane
- [ ] Zadanie do zrobienia

### Tabela
| Cecha | Wartość | Status |
| :--- | :--- | :--- |
| Test 1 | Ok | Zaliczone |
| Test 2 | Błąd | Do poprawy |

### Kod w języku Python
```python
import numpy as np

def analizuj_wektory(p1_poczatek, p1_koniec, p2_poczatek, p2_koniec):
    """
    Oblicza długości, iloczyn skalarny oraz iloczyn wektorowy 
    dla dwóch wektorów zdefiniowanych punktami początkowymi i końcowymi (2D lub 3D).
    """
    # Konwersja punktów na tablice numpy
    A_pocz = np.array(p1_poczatek, dtype=float)
    A_kon = np.array(p1_koniec, dtype=float)
    B_pocz = np.array(p2_poczatek, dtype=float)
    B_kon = np.array(p2_koniec, dtype=float)
    
    # Wyznaczenie składowych wektorów: v = K - P
    v1 = A_kon - A_pocz
    v2 = B_kon - B_pocz
    
    # Obliczenie długości wektorów (norma Euklidesowa)
    dlugosc_v1 = np.linalg.norm(v1)
    dlugosc_v2 = np.linalg.norm(v2)
    
    # Iloczyn skalarny
    iloczyn_skalarny = np.dot(v1, v2)
    
    # Iloczyn wektorowy (dla 2D uzupełniamy z-ową współrzędną zerem)
    v1_3d = np.pad(v1, (0, 3 - len(v1))) if len(v1) == 2 else v1
    v2_3d = np.pad(v2, (0, 3 - len(v2))) if len(v2) == 2 else v2
    
    iloczyn_wektorowy = np.cross(v1_3d, v2_3d)
    
    return {
        "wektor_1": v1,
        "dlugosc_1": dlugosc_v1,
        "wektor_2": v2,
        "dlugosc_2": dlugosc_v2,
        "iloczyn_skalarny": iloczyn_skalarny,
        "iloczyn_wektorowy_3D": iloczyn_wektorowy
    }


# === PRZYKŁAD UŻYCIA 2D ===
print("--- PRZYKŁAD 2D ---")
p1_start = [1, 2]
p1_end   = [4, 6]
p2_start = [0, 0]
p2_end   = [3, 1]

wyniki_2d = analizuj_wektory(p1_start, p1_end, p2_start, p2_end)
print(f"Wektor 1: {wyniki_2d['wektor_1']}, Długość: {wyniki_2d['dlugosc_1']:.2f}")
print(f"Wektor 2: {wyniki_2d['wektor_2']}, Długość: {wyniki_2d['dlugosc_2']:.2f}")
print(f"Iloczyn skalarny: {wyniki_2d['iloczyn_skalarny']}")
print(f"Iloczyn wektorowy (wynik w 3D): {wyniki_2d['iloczyn_wektorowy_3D']}\n")


# === PRZYKŁAD UŻYCIA 3D ===
print("--- PRZYKŁAD 3D ---")
A_pocz = [1, 0, 2]
A_kon  = [3, 1, 5]
B_pocz = [2, -1, 1]
B_kon  = [0, 4, 2]

wyniki_3d = analizuj_wektory(A_pocz, A_kon, B_pocz, B_kon)
print(f"Wektor 1: {wyniki_3d['wektor_1']}, Długość: {wyniki_3d['dlugosc_1']:.2f}")
print(f"Wektor 2: {wyniki_3d['wektor_2']}, Długość: {wyniki_3d['dlugosc_2']:.2f}")
print(f"Iloczyn skalarny: {wyniki_3d['iloczyn_skalarny']}")
print(f"Iloczyn wektorowy: {wyniki_3d['iloczyn_wektorowy_3D']}")
