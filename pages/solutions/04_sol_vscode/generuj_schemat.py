from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

# Konfiguracja ścieżek i stylu
OUT = Path("files/rzut_ukosny")
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, INK = "#2463A6", "#D47722", "#233044"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.labelcolor": INK,
    "text.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def zapisz(fig, nazwa):
    """Zapisuje wykres w wysokiej rozdzielczości i wyświetla go."""
    fig.savefig(OUT / f"{nazwa}.png", dpi=300, bbox_inches="tight")
    if plt.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)


# --- Dane do schematu (poglądowe) ---
v0_mag = 5.0
alpha_deg = 45.0
alpha_rad = np.radians(alpha_deg)
g = 9.81

# Całkowity czas lotu i punkt pokazowy na torze (35% czasu lotu)
t_calkowity = (2 * v0_mag * np.sin(alpha_rad)) / g
czas_pokazu = 0.35 * t_calkowity

x_pokazu = v0_mag * np.cos(alpha_rad) * czas_pokazu
y_pokazu = (v0_mag * np.sin(alpha_rad) * czas_pokazu) - (0.5 * g * czas_pokazu**2)

# Składowe prędkości w punkcie pokazowym
vx_pokazu = v0_mag * np.cos(alpha_rad)
vy_pokazu = (v0_mag * np.sin(alpha_rad)) - (g * czas_pokazu)

# --- Generowanie rysunku ---
fig, ax = plt.subplots(figsize=(10, 6))

# 1. Ustawienia osi
ax.set(xlim=(-0.6, v0_mag + 0.8), ylim=(-0.6, y_pokazu + 2.2))
ax.set_aspect("equal")
ax.set_xlabel("Oś pozioma x [m]", loc="right")
ax.set_ylabel("Oś pionowa y [m]", loc="top")

ax.spines["left"].set_position("zero")
ax.spines["bottom"].set_position("zero")
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)
ax.text(-0.2, -0.25, "0", fontsize=12, ha="right")

# 2. Tor lotu (linia przerywana)
t_lin = np.linspace(0, t_calkowity, 500)
x_lin = v0_mag * np.cos(alpha_rad) * t_lin
y_lin = (v0_mag * np.sin(alpha_rad) * t_lin) - (0.5 * g * t_lin**2)
ax.plot(x_lin, y_lin, color=BLUE, linestyle="--", alpha=0.7)

# 3. Wektor prędkości początkowej v0
ax.annotate(
    "",
    xy=(v0_mag * np.cos(alpha_rad) * 0.45, v0_mag * np.sin(alpha_rad) * 0.45),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2.5),
)
ax.text(
    v0_mag * np.cos(alpha_rad) * 0.22 - 0.15,
    v0_mag * np.sin(alpha_rad) * 0.22 + 0.2,
    r"$\vec{v}_0$",
    color=ORANGE,
    fontsize=14,
    fontweight="bold",
)

# 4. Łuk kąta alfa
ax.add_patch(Arc((0, 0), 1.4, 1.4, theta1=0, theta2=alpha_deg, color=INK, lw=1.2))
ax.text(0.85, 0.3, r"$\alpha$", fontsize=12, color=INK)

# 5. Ciało w locie (kropka)
ax.plot(x_pokazu, y_pokazu, "o", color=ORANGE, markersize=10, zorder=5)

# 6. Wektory sił i prędkości w punkcie pokazowym
# Prędkość v (styczna do toru, pomarańczowa)
ax.arrow(
    x_pokazu,
    y_pokazu,
    vx_pokazu * 0.28,
    vy_pokazu * 0.28,
    head_width=0.12,
    head_length=0.16,
    fc=ORANGE,
    ec=ORANGE,
    lw=2,
    zorder=4,
)
ax.text(
    x_pokazu + vx_pokazu * 0.28 + 0.15,
    y_pokazu + vy_pokazu * 0.28 + 0.05,
    r"$\vec{v}$",
    color=ORANGE,
    fontsize=14,
    fontweight="bold",
)

# Grawitacja Fg (pionowo w dół, niebieska)
dlugosc_wektora_g = 1.1
ax.arrow(
    x_pokazu,
    y_pokazu,
    0,
    -dlugosc_wektora_g,
    head_width=0.12,
    head_length=0.16,
    fc=BLUE,
    ec=BLUE,
    lw=2,
    zorder=4,
)
ax.text(
    x_pokazu + 0.15,
    y_pokazu - dlugosc_wektora_g * 0.6,
    r"$\vec{F}_g = m\vec{g}$",
    color=BLUE,
    fontsize=12,
)

# Opór powietrza Fop (przeciwny do v, niebieski)
ax.arrow(
    x_pokazu,
    y_pokazu,
    -vx_pokazu * 0.22,
    -vy_pokazu * 0.22,
    head_width=0.12,
    head_length=0.16,
    fc=BLUE,
    ec=BLUE,
    lw=2,
    zorder=4,
)
ax.text(
    x_pokazu - vx_pokazu * 0.22 - 0.55,
    y_pokazu - vy_pokazu * 0.22 - 0.15,
    r"$\vec{F}_{\mathrm{op}}$",
    color=BLUE,
    fontsize=12,
)

# Zapisanie obrazu
zapisz(fig, "schemat")
print("Schemat wygenerowany poprawnie przy użyciu Pythona.")