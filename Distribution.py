from _generator import generate
import matplotlib.pyplot as plt

# Generar una muestra de control con señal J/psi y fondo
df = generate(
    sample="control",
    yields={"bkg": 10000, "jpsi": 10000},
    seed=42
)

# Resumen estadístico de las variables
columns = ["mu1_pid", "mu2_pid", "vtx_chi2"]

print(df[columns].describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9]))

# Crear una figura con 1 fila y 3 columnas
fig, ax = plt.subplots(1, 3, figsize=(15, 5))

# Distribución del PID del primer muón
ax[0].hist(df["mu1_pid"].dropna(), bins=50)
ax[0].set_xlabel("mu1_pid")
ax[0].set_ylabel("Número de candidatos")
ax[0].set_title("Identificación del muón 1")

# Distribución del PID del segundo muón
ax[1].hist(df["mu2_pid"].dropna(), bins=50)
ax[1].set_xlabel("mu2_pid")
ax[1].set_ylabel("Número de candidatos")
ax[1].set_title("Identificación del muón 2")

# Distribución de la calidad del vértice
ax[2].hist(df["vtx_chi2"].dropna(), bins=50)
ax[2].set_xlabel("vtx_chi2")
ax[2].set_ylabel("Número de eventos")
ax[2].set_title("Calidad del vértice")

# Ajustar espacios y mostrar la figura
plt.tight_layout()
plt.show()