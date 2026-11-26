from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from scipy.integrate import solve_ivp

# =====================================================================
# PARAMETRY Z NUMERU INDEKSU (UZUPEŁNIJ PRZED URUCHOMIENIEM)
# =====================================================================
V0 = 36.0              
ALPHA_DEG = 70.0       
POKAZ_WYKRESY = False  # Zostaw False, żeby zapisać wszystkie naraz bez zatrzymywania

if V0 <= 0:
    raise ValueError("Prędkość początkowa V0 musi być dodatnia.")
if not (0.0 < ALPHA_DEG < 90.0):
    raise ValueError("Kąt wyrzutu ALPHA_DEG musi należeć do przedziału (0, 90) stopni.")

# Parametry fizyczne w jednostkach SI
g = 9.81                       # przyspieszenie ziemskie [m/s^2]
m = 0.20                       # masa ciała [kg]
R = 0.035                      # promień kuli [m]
rho = 1.225                    # gęstość powietrza [kg/m^3]
Cd = 0.47                      # współczynnik oporu kuli [-]
S = np.pi * R**2               # pole przekroju [m^2]
c = 0.5 * rho * Cd * S         # współczynnik oporu kwadratowego [kg/m]
b = 0.050                      # współczynnik oporu liniowego [kg/s]
gamma = b / m                  # parametr tłumienia liniowego [1/s]

alpha = np.deg2rad(ALPHA_DEG)
v0x = V0 * np.cos(alpha)
v0y = V0 * np.sin(alpha)

BASE = Path(__file__).resolve().parent if "__file__" in globals() else Path(".")
OUT = BASE / "files" / "rzut_ukosny"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, GREEN = "#2463A6", "#D47722", "#168578"
PURPLE, INK = "#8455A4", "#233044"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.labelcolor": INK,
    "text.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.18,
    "lines.linewidth": 2.3,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def zapisz(fig, nazwa, dpi=300):
    """Zapisuje wykres w wysokiej rozdzielczości i opcjonalnie go wyświetla."""
    sciezka = OUT / f"{nazwa}.png"
    fig.savefig(sciezka, dpi=dpi, bbox_inches="tight")
    print(f"Zapisano: {sciezka.name}")
    if POKAZ_WYKRESY and plt.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------
# BLOK A: Schemat układu fizycznego (Rysunek 1)
# ---------------------------------------------------------------------
tw_id = v0y / g
tc_id = 2.0 * tw_id
H_id = v0y**2 / (2.0 * g)
Z_id = V0**2 * np.sin(2.0 * alpha) / g
E0 = 0.5 * m * V0**2

v0_mag = 5.0
alpha_sch_deg = 45.0
alpha_sch_rad = np.radians(alpha_sch_deg)

t_calkowity = (2 * v0_mag * np.sin(alpha_sch_rad)) / g
czas_pokazu = 0.35 * t_calkowity

x_pokazu = v0_mag * np.cos(alpha_sch_rad) * czas_pokazu
y_pokazu = (v0_mag * np.sin(alpha_sch_rad) * czas_pokazu) - (0.5 * g * czas_pokazu**2)

vx_pokazu = v0_mag * np.cos(alpha_sch_rad)
vy_pokazu = (v0_mag * np.sin(alpha_sch_rad)) - (g * czas_pokazu)

fig, ax = plt.subplots(figsize=(10, 6))
ax.grid(False)

ax.set(xlim=(-0.6, v0_mag + 0.8), ylim=(-0.6, y_pokazu + 2.2))
ax.set_aspect("equal")
ax.set_xlabel("Oś pozioma x [m]", loc="right", fontsize=12)
ax.set_ylabel("Oś pionowa y [m]", loc="top", fontsize=12)

ax.spines["left"].set_position("zero")
ax.spines["bottom"].set_position("zero")
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)
ax.text(-0.2, -0.25, "0", fontsize=12, ha="right")

t_sch = np.linspace(0, t_calkowity, 500)
x_sch = v0_mag * np.cos(alpha_sch_rad) * t_sch
y_sch = (v0_mag * np.sin(alpha_sch_rad) * t_sch) - (0.5 * g * t_sch**2)
ax.plot(x_sch, y_sch, color=BLUE, linestyle="--", alpha=0.7)

ax.annotate(
    "",
    xy=(v0_mag * np.cos(alpha_sch_rad) * 0.45, v0_mag * np.sin(alpha_sch_rad) * 0.45),
    xytext=(0, 0),
    arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2.5),
)
ax.text(
    v0_mag * np.cos(alpha_sch_rad) * 0.22 - 0.15,
    v0_mag * np.sin(alpha_sch_rad) * 0.22 + 0.2,
    r"$\vec{v}_0$",
    color=ORANGE,
    fontsize=14,
    fontweight="bold",
)

ax.add_patch(Arc((0, 0), 1.4, 1.4, theta1=0, theta2=alpha_sch_deg, color=INK, lw=1.2))
ax.text(0.85, 0.3, r"$\alpha$", fontsize=12, color=INK)

ax.plot(x_pokazu, y_pokazu, "o", color=ORANGE, markersize=10, zorder=5)

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

zapisz(fig, "schemat")

# ---------------------------------------------------------------------
# BLOK B: Ruch idealny i energia (Rysunki 2 i 3)
# ---------------------------------------------------------------------
t_id = np.linspace(0.0, tc_id, 1001)
x_id = v0x * t_id
y_id = v0y * t_id - 0.5 * g * t_id**2
vx_id = np.full_like(t_id, v0x)
vy_id = v0y - g * t_id
v_id = np.hypot(vx_id, vy_id)

Ep_id = m * g * y_id
Ek_id = 0.5 * m * v_id**2

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
axes[0].plot(x_id, y_id, color=BLUE, label="Tor $y(x)$")
axes[0].scatter([Z_id / 2, Z_id], [H_id, 0.0], color=ORANGE, zorder=4)
axes[0].axvline(Z_id / 2, color=INK, ls=":", lw=1.0, alpha=0.5)
axes[0].annotate(
    f"$H_{{max}} = {H_id:.2f}$ m",
    xy=(Z_id / 2, H_id),
    xytext=(Z_id / 2 + 0.06 * Z_id, H_id * 0.92),
    arrowprops={"arrowstyle": "->", "color": INK},
)
axes[0].annotate(
    f"$Z = {Z_id:.2f}$ m",
    xy=(Z_id, 0),
    xytext=(Z_id * 0.72, H_id * 0.18),
    arrowprops={"arrowstyle": "->", "color": INK},
)
axes[0].set(
    xlabel="Położenie poziome x [m]",
    ylabel="Wysokość y [m]",
    title="Tor lotu bez oporu powietrza",
    ylim=(0, H_id * 1.15),
)
axes[0].legend(loc="upper right")

axes[1].plot(t_id, vx_id, color=BLUE, label="$v_x(t)$")
axes[1].plot(t_id, vy_id, color=ORANGE, label="$v_y(t)$")
axes[1].plot(t_id, v_id, color=GREEN, ls="--", label="$v(t)$")
axes[1].axhline(0, color=INK, lw=0.8, alpha=0.5)
axes[1].axvline(tw_id, color=INK, ls=":", lw=1.0, alpha=0.5, label="$t_w$")
axes[1].set(
    xlabel="Czas t [s]",
    ylabel="Prędkość [m/s]",
    title="Składowe prędkości i szybkość",
)
axes[1].legend(loc="lower left", ncol=2, fontsize=9)
fig.tight_layout()
zapisz(fig, "tor_i_predkosc_idealny")

fig, ax = plt.subplots(figsize=(10, 4.2))
ax.plot(t_id, Ek_id, color=ORANGE, label="Kinetyczna $E_k$")
ax.plot(t_id, Ep_id, color=BLUE, label="Potencjalna $E_p$")
ax.plot(t_id, Ek_id + Ep_id, color=GREEN, ls="--", label="Całkowita $E$")
ax.axvline(tw_id, color=INK, ls=":", lw=1.0, alpha=0.5)
ax.set(xlabel="Czas t [s]", ylabel="Energia [J]", ylim=(-0.03 * E0, 1.15 * E0))
ax.set_title("Zachowanie energii mechanicznej w rzucie bez oporu", loc="left")
ax.legend(loc="center right", fontsize=9)
fig.tight_layout()
zapisz(fig, "energia_idealna")

assert np.max(np.abs(Ek_id + Ep_id - E0)) < 1e-11


# ---------------------------------------------------------------------
# BLOK C: Opór liniowy i kwadratowy (Rysunki 4 i 5)
# ---------------------------------------------------------------------
def upadek(t, stan):
    return stan[1]


upadek.terminal = True
upadek.direction = -1


def symuluj_rzut(v0_val, alpha_rad, tryb="kwadratowy", b_val=b, c_val=c, m_val=m):
    vx0 = v0_val * np.cos(alpha_rad)
    vy0 = v0_val * np.sin(alpha_rad)
    t_max = max(6.0 * (2.0 * vy0 / g) + 30.0, 60.0)

    def rhs(t, stan):
        xx, yy, vx, vy = stan
        if tryb == "brak":
            ax_val, ay_val = 0.0, -g
        elif tryb == "liniowy":
            ax_val = -(b_val / m_val) * vx
            ay_val = -g - (b_val / m_val) * vy
        elif tryb == "kwadratowy":
            vv = np.hypot(vx, vy)
            ax_val = -(c_val / m_val) * vv * vx
            ay_val = -g - (c_val / m_val) * vv * vy
        else:
            raise ValueError("Nieznany tryb oporu.")
        return [vx, vy, ax_val, ay_val]

    sol = solve_ivp(
        rhs,
        (0.0, t_max),
        [0.0, 0.0, vx0, vy0],
        events=upadek,
        dense_output=True,
        rtol=1e-9,
        atol=1e-11,
    )
    if not sol.success or len(sol.t_events[0]) == 0:
        raise RuntimeError("Nie udało się wyznaczyć chwili upadku.")

    t_hit = sol.t_events[0][0]
    ts = np.linspace(0.0, t_hit, 1001)
    xs_s, ys_s, vxs_s, vys_s = sol.sol(ts)
    ys_s = np.maximum(ys_s, 0.0)
    return ts, xs_s, ys_s, vxs_s, vys_s


t_lin, x_lin, y_lin, vx_lin, vy_lin = symuluj_rzut(V0, alpha, "liniowy")
t_kw, x_kw, y_kw, vx_kw, vy_kw = symuluj_rzut(V0, alpha, "kwadratowy")

v_gr = g / gamma
x_lin_an = (v0x / gamma) * (1.0 - np.exp(-gamma * t_lin))
y_lin_an = ((v0y + v_gr) / gamma) * (1.0 - np.exp(-gamma * t_lin)) - v_gr * t_lin
assert np.max(np.abs(x_lin - x_lin_an)) < 1e-7
assert np.max(np.abs(y_lin - y_lin_an)) < 1e-7

v_lin = np.hypot(vx_lin, vy_lin)
v_kw = np.hypot(vx_kw, vy_kw)
E_lin = 0.5 * m * v_lin**2 + m * g * y_lin
E_kw = 0.5 * m * v_kw**2 + m * g * y_kw
assert np.max(np.diff(E_lin)) < 1e-10 and np.max(np.diff(E_kw)) < 1e-10

fig, ax = plt.subplots(figsize=(10, 4.6))
ax.plot(x_id, y_id, "--", color=BLUE, label="Próżnia (model idealny)")
ax.plot(x_lin, y_lin, color=GREEN, label=f"Opór liniowy ($b = {b:.3f}$ kg/s)")
ax.plot(x_kw, y_kw, color=ORANGE, label=f"Opór kwadratowy ($c = {c:.4f}$ kg/m)")
ax.set(
    xlabel="Odległość pozioma x [m]",
    ylabel="Wysokość y [m]",
    xlim=(0, Z_id * 1.05),
    ylim=(0, H_id * 1.15),
)
ax.set_title("Porównanie trajektorii lotu w trzech modelach", loc="left")
ax.legend(loc="upper right")
fig.tight_layout()
zapisz(fig, "porownanie_torow")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
axes[0].plot(t_id, Ek_id + Ep_id, "--", color=BLUE, label="Próżnia")
axes[0].plot(t_lin, E_lin, color=GREEN, label="Opór liniowy")
axes[0].plot(t_kw, E_kw, color=ORANGE, label="Opór kwadratowy")
axes[0].set(
    xlabel="Czas t [s]",
    ylabel="Energia całkowita E [J]",
    title="Spadek energii mechanicznej",
    ylim=(0, E0 * 1.12),
)
axes[0].legend(loc="lower left")

axes[1].plot(t_id, v_id, "--", color=BLUE, label="Próżnia")
axes[1].plot(t_lin, v_lin, color=GREEN, label="Opór liniowy")
axes[1].plot(t_kw, v_kw, color=ORANGE, label="Opór kwadratowy")
axes[1].set(
    xlabel="Czas t [s]",
    ylabel="Szybkość v [m/s]",
    title="Zmiana szybkości podczas lotu",
    ylim=(0, V0 * 1.08),
)
axes[1].legend(loc="lower left")
fig.tight_layout()
zapisz(fig, "energia_i_predkosc_opor")

# ---------------------------------------------------------------------
# BLOK D: Optymalny kąt rzutu i wpływ oporu (Rysunki 6 i 7)
# ---------------------------------------------------------------------
katy_deg = np.linspace(5.0, 85.0, 81)
katy_rad = np.deg2rad(katy_deg)

fig, ax = plt.subplots(figsize=(10, 4.5))
Z_proznia = (V0**2 * np.sin(2.0 * katy_rad)) / g
ax.plot(katy_deg, Z_proznia, "--", color=BLUE, label="Próżnia (opt. 45.0°)")

for mnoznik, color in [(0.5, GREEN), (1.0, ORANGE), (3.0, PURPLE)]:
    c_test = mnoznik * c
    zasiegi = np.array([
        symuluj_rzut(V0, k_rad, "kwadratowy", c_val=c_test)[1][-1]
        for k_rad in katy_rad
    ])
    idx_max = int(np.argmax(zasiegi))
    kat_opt = float(katy_deg[idx_max])
    ax.plot(
        katy_deg,
        zasiegi,
        color=color,
        label=f"$c = {c_test:.4f}$ kg/m (opt. {kat_opt:.1f}°)",
    )
    ax.scatter([kat_opt], [zasiegi[idx_max]], color=color, zorder=5)

ax.axvline(45.0, color=INK, ls=":", lw=1.0, alpha=0.5)
ax.set(xlabel="Kąt wyrzutu [°]", ylabel="Zasięg Z [m]", xlim=(5, 85))
ax.set_title("Zależność zasięgu od kąta wyrzutu — przesunięcie optimum", loc="left")
ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
zapisz(fig, "optymalny_kat")

fig, ax = plt.subplots(figsize=(10, 4.5))
ax.plot(x_id, y_id, "--", color=BLUE, label="$c/m = 0$ (próżnia)")
for mnoznik, color in [(0.5, GREEN), (1.0, ORANGE), (3.0, PURPLE), (8.0, INK)]:
    c_test = mnoznik * c
    _, xx, yy, _, _ = symuluj_rzut(V0, alpha, "kwadratowy", c_val=c_test)
    ax.plot(xx, yy, color=color, label=f"$c/m = {c_test/m:.4f}$ m⁻¹")

ax.set(
    xlabel="Odległość pozioma x [m]",
    ylabel="Wysokość y [m]",
    xlim=(0, Z_id * 1.05),
    ylim=(0, H_id * 1.15),
)
ax.set_title("Asymetria krzywej balistycznej wraz ze wzrostem oporu c/m", loc="left")
ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
zapisz(fig, "wplyw_oporu")

print(f"Gotowe: zapisano 7 rysunków w {OUT.resolve()}")