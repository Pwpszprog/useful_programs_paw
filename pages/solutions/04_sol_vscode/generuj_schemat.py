from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from scipy.integrate import solve_ivp

# =====================================================================
# PARAMETRY Z NUMERU INDEKSU (UZUPEŁNIJ PRZED URUCHOMIENIEM)
# =====================================================================
V0 = 36.0              # <-- Dwie pierwsze cyfry indeksu [m/s]
ALPHA_DEG = 70.0       # <-- Dwie ostatnie cyfry indeksu [stopnie]
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

alpha = np.deg2rad(ALPHA_DEG)  # kąt w jednostkach SI [rad]
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


def zapisz(fig, nazwa, dpi=180):
    """Zapisuje wykres w rozdzielczości 180 dpi i opcjonalnie go wyświetla."""
    sciezka = OUT / f"{nazwa}.png"
    fig.savefig(sciezka, dpi=dpi, bbox_inches="tight")
    print(f"Zapisano: {sciezka.name}")
    if POKAZ_WYKRESY and plt.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------
# BLOK A: Parametry modelu idealnego i schemat układu fizycznego (Rysunek 1)
# ---------------------------------------------------------------------
tw_id = v0y / g
tc_id = 2.0 * tw_id
H_id = v0y**2 / (2.0 * g)
Z_id = V0**2 * np.sin(2.0 * alpha) / g
E0 = 0.5 * m * V0**2

print("=== WYNIKI DO ROZDZIAŁU 1 i 3 (MODEL IDEALNY, SI) ===")
print(f"v0 = {V0:.2f} m/s | alpha = {ALPHA_DEG:.1f} deg = {alpha:.4f} rad")
print(f"v0x = {v0x:.3f} m/s | v0y = {v0y:.3f} m/s")
print(f"tw = {tw_id:.3f} s | tc = {tc_id:.3f} s")
print(f"H_max = {H_id:.3f} m | Z = {Z_id:.3f} m | E = {E0:.3f} J")
print(f"c = {c:.5f} kg/m | v_gr_kw = {np.sqrt(m * g / c):.2f} m/s\n")

fig, ax = plt.subplots(figsize=(10, 4.2), layout="constrained")
ax.set(xlim=(-0.6, 8.2), ylim=(-0.6, 3.6))
ax.axis("off")

# Osie układu współrzędnych
ax.annotate("", xy=(7.8, 0), xytext=(-0.3, 0),
            arrowprops={"arrowstyle": "->", "color": INK, "lw": 1.6})
ax.annotate("", xy=(0, 3.3), xytext=(0, -0.3),
            arrowprops={"arrowstyle": "->", "color": INK, "lw": 1.6})
ax.text(7.85, -0.1, "$x$", fontsize=13)
ax.text(-0.25, 3.25, "$y$", fontsize=13)
ax.text(-0.35, -0.3, "$(0,0)$", fontsize=10)

# Przykładowe tory (schematycznie)
xs = np.linspace(0, 7.0, 300)
ys_id = 1.45 * xs * (1 - xs / 7.0)
ys_op = 1.45 * xs * (1 - xs / 5.4) * np.exp(-0.10 * xs)
mask = ys_op >= 0
ax.plot(xs, ys_id, "--", color=BLUE, alpha=0.75, label="Tor idealny (parabola)")
ax.plot(xs[mask], ys_op[mask], color=ORANGE, label="Tor z oporem powietrza")

# Wektor prędkości początkowej i kąt alfa
ax.annotate("", xy=(1.5, 2.18), xytext=(0, 0),
            arrowprops={"arrowstyle": "->", "color": GREEN, "lw": 2.5})
ax.text(1.55, 2.25, r"$\vec{v}_0$", color=GREEN, fontsize=14, fontweight="bold")
ax.add_patch(Arc((0, 0), 1.8, 1.8, theta1=0, theta2=np.rad2deg(np.arctan(1.45)),
                 color=INK, lw=1.4))
ax.text(1.02, 0.42, r"$\alpha$", fontsize=13)

# Punkt na torze z wektorami sił
idx = 70
px, py = xs[idx], ys_op[idx]
tx = xs[idx + 1] - xs[idx - 1]
ty = ys_op[idx + 1] - ys_op[idx - 1]
tnorm = np.hypot(tx, ty)
ux, uy = tx / tnorm, ty / tnorm

ax.scatter([px], [py], color=INK, s=75, zorder=5)
ax.annotate("", xy=(px + 1.15 * ux, py + 1.15 * uy), xytext=(px, py),
            arrowprops={"arrowstyle": "->", "color": GREEN, "lw": 2.2})
ax.text(px + 1.2 * ux, py + 1.2 * uy + 0.08, r"$\vec{v}$", color=GREEN, fontsize=13)

ax.annotate("", xy=(px, py - 1.15), xytext=(px, py),
            arrowprops={"arrowstyle": "->", "color": BLUE, "lw": 2.2})
ax.text(px + 0.12, py - 1.1, r"$\vec{F}_g = m\vec{g}$", color=BLUE, fontsize=12)

ax.annotate("", xy=(px - 0.95 * ux, py - 0.95 * uy), xytext=(px, py),
            arrowprops={"arrowstyle": "->", "color": PURPLE, "lw": 2.2})
ax.text(px - 1.05 * ux - 0.45, py - 0.95 * uy - 0.2,
        r"$\vec{F}_{\mathrm{op}}$", color=PURPLE, fontsize=12)

ax.set_title("Model: rzut ukośny w polu grawitacyjnym (próżnia i opór powietrza)", loc="left")
ax.legend(loc="upper right", frameon=True)
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

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout="constrained")
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
zapisz(fig, "tor_i_predkosc_idealny")

fig, ax = plt.subplots(figsize=(10, 4.2), layout="constrained")
ax.plot(t_id, Ek_id, color=ORANGE, label="Kinetyczna $E_k$")
ax.plot(t_id, Ep_id, color=BLUE, label="Potencjalna $E_p$")
ax.plot(t_id, Ek_id + Ep_id, color=GREEN, ls="--", label="Całkowita $E$")
ax.axvline(tw_id, color=INK, ls=":", lw=1.0, alpha=0.5)
ax.set(xlabel="Czas t [s]", ylabel="Energia [J]", ylim=(-0.03 * E0, 1.15 * E0))
ax.set_title("Zachowanie energii mechanicznej w rzucie bez oporu", loc="left")
ax.legend(loc="center right", fontsize=9)
zapisz(fig, "energia_idealna")

blad_E = np.max(np.abs(Ek_id + Ep_id - E0))
assert blad_E < 1e-11
print(f"Kontrola energii (model idealny): max błąd = {blad_E:.3e} J\n")


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
blad_lin_x = np.max(np.abs(x_lin - x_lin_an))
blad_lin_y = np.max(np.abs(y_lin - y_lin_an))
assert blad_lin_x < 1e-7 and blad_lin_y < 1e-7
print(f"Kontrola (opór liniowy): błąd x = {blad_lin_x:.3e} m, błąd y = {blad_lin_y:.3e} m")

v_lin = np.hypot(vx_lin, vy_lin)
v_kw = np.hypot(vx_kw, vy_kw)
E_lin = 0.5 * m * v_lin**2 + m * g * y_lin
E_kw = 0.5 * m * v_kw**2 + m * g * y_kw
assert np.max(np.diff(E_lin)) < 1e-10 and np.max(np.diff(E_kw)) < 1e-10
print("Kontrola dysypacji: energia w obu modelach z oporem nie rośnie.\n")

print("=== WYNIKI DO TABELI W ROZDZIALE 6.1 ===")
for nazwa, tt, xx, yy, vv, EE in [
    ("Próżnia", t_id, x_id, y_id, v_id, Ek_id + Ep_id),
    ("Opór liniowy", t_lin, x_lin, y_lin, v_lin, E_lin),
    ("Opór kwadratowy", t_kw, x_kw, y_kw, v_kw, E_kw),
]:
    strata = (1.0 - EE[-1] / E0) * 100.0
    print(f"{nazwa:16s} | Z = {xx[-1]:6.2f} m | H_max = {np.max(yy):6.2f} m | "
          f"tc = {tt[-1]:5.2f} s | v(tc) = {vv[-1]:5.2f} m/s | strata E = {strata:5.1f} %")
print()

fig, ax = plt.subplots(figsize=(10, 4.6), layout="constrained")
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
zapisz(fig, "porownanie_torow")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout="constrained")
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
zapisz(fig, "energia_i_predkosc_opor")

# ---------------------------------------------------------------------
# BLOK D: Optymalny kąt rzutu i wpływ oporu (Rysunki 6 i 7)
# ---------------------------------------------------------------------
katy_deg = np.linspace(5.0, 85.0, 81)
katy_rad = np.deg2rad(katy_deg)

fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
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
zapisz(fig, "optymalny_kat")

fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
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
zapisz(fig, "wplyw_oporu")

print(f"Gotowe: zapisano 7 rysunków w {OUT.resolve()}")
