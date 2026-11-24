import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Konfiguracja ścieżek i stylu
OUT = Path("files/rzut_ukosny")
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, INK = "#2463A6", "#D47722", "#233044"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
})

defzapisz(fig, nazwa):
    """Zapisuje wykres w wysokiej rozdzielczości i wyświetla go."""
    fig.savefig(OUT / f"{nazwa}.png", dpi=300, bbox_inches='tight')
    if plt.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)

# --- Dane do schematu (poglądowe) ---
v0_mag = 5.0
alpha_deg = 45.0
alpha_rad = np.radians(alpha_deg)
g = 9.81

# Obliczenie punktu w połowie drogi (dla wizualizacji wektorów)
# Zakładamy lot do upadku na y=0 w modelu idealnym (opór tylko dorysowujemy)
t_polowa = (2 * v0_mag * np.sin(alpha_rad)) / g
czas_pokazu = t_polowa / 2 

x_pokazu = v0_mag * np.cos(alpha_rad) * czas_pokazu
y_pokazu = (v0_mag * np.sin(alpha_rad) * czas_pokazu) - (0.5 * g * czas_pokazu**2)

# Składowe prędkości w tym punkcie (idealnej)
vx_pokazu = v0_mag * np.cos(alpha_rad)
vy_pokazu = (v0_mag * np.sin(alpha_rad)) - (g * czas_pokazu)
v_pokazu = np.hypot(vx_pokazu, vy_pokazu)

# --- Generowanie rysunku ---
fig, ax = plt.subplots(figsize=(10, 6))

# 1. Ustawienia osi
ax.set(xlim=(-1, v0_mag + 1), ylim=(-1, y_pokazu + 2))
ax.set_aspect('equal')
ax.set_xlabel("Oś pozioma $x$ [m]")
ax.set_ylabel("Oś pionowa $y$ [m]")
# Ukrywamy domyślne osie, dorysowujemy własne strzałki
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)
ax.text(-0.5, 0, "0", fontsize=12, ha='right') # Etykieta początku układu

# 2. Tor lotu (linia przerywana)
# Generujemy gładką parabolę (trajektoria idealna)
t_max = (2 * v0_mag * np.sin(alpha_rad)) / g
t_lin = np.linspace(0, t_max, 500)
x_lin = v0_mag * np.cos(alpha_rad) * t_lin
y_lin = (v0_mag * np.sin(alpha_rad) * t_lin) - (0.5 * g * t_lin**2)
ax.plot(x_lin, y_lin, color=BLUE, linestyle='--', alpha=0.7, label='Tor lotu (trajektoria)')

# 3. Wektor prędkości początkowej v0
ax.annotate('', xy=(v0_mag * np.cos(alpha_rad), v0_mag * np.sin(alpha_rad)), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color=ORANGE, lw=2.5))
ax.text(v0_mag * np.cos(alpha_rad) * 0.5 + 0.2, v0_mag * np.sin(alpha_rad) * 0.5 + 0.2, 
        r'$\vec{v}_0$', color=ORANGE, fontsize=14, fontweight='bold')

# 4. Łuk kąta alfa
from matplotlib.patches import Arc
ax.add_patch(Arc((0, 0), 2, 2, theta1=0, theta2=alpha_deg, color=INK, lw=1.2))
ax.text(1.2, 0.5, r'$\alpha$', fontsize=12, color=INK)

# 5. Ciało w locie (kropka)
ax.plot(x_pokazu, y_pokazu, 'o', color=ORANGE, markersize=10)

# 6. Wektory sił i prędkości w punkcie pokazowym
# Prędkość v (styczna, pomarańczowa)
ax.arrow(x_pokazu, y_pokazu, vx_pokazu*0.3, vy_pokazu*0.3, 
         head_width=0.15, head_length=0.2, fc=ORANGE, ec=ORANGE, lw=2)
ax.text(x_pokazu + vx_pokazu*0.3 + 0.2, y_pokazu + vy_pokazu*0.3, 
        r'$\vec{v}$', color=ORANGE, fontsize=14, fontweight='bold')

# Grawitacja Fg (w dół, niebieska)
dlugosc_wektora_g = 1.5
ax.arrow(x_pokazu, y_pokazu, 0, -dlugosc_wektora_g, 
         head_width=0.15, head_length=0.2, fc=BLUE, ec=BLUE, lw=2)
ax.text(x_pokazu + 0.2, y_pokazu - dlugosc_wektora_g / 2, 
        r'$\vec{F}_g = m\vec{g}$', color=BLUE, fontsize=12)

# Opór powietrza Fop (przeciwny do v, niebieski)
ax.arrow(x_pokazu, y_pokazu, -vx_pokazu*0.2, -vy_pokazu*0.2, 
         head_width=0.15, head_length=0.2, fc=BLUE, ec=BLUE, lw=2)
ax.text(x_pokazu - vx_pokazu*0.2 - 0.8, y_pokazu - vy_pokazu*0.2, 
        r'$\vec{F}_{\mathrm{op}}$', color=BLUE, fontsize=12)

# 7. Tytuł i legenda (opcjonalne, w raporcie podpisujemy rysunek)
# ax.set_title("Schemat sił działających w rzucie ukośnym", fontweight='bold')
# ax.legend(loc='upper right', frameon=False)

# Zapisanie obrazu
zapisz(fig, "schemat")
print("Schemat wygenerowany poprawnie przy użyciu Pythona.")