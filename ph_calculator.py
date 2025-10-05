
# ph_calculator.py

import math

# === Константы ===
MOLAR_MASSES = {
    'HCOONa': 68.01,   # формиат натрия
    'HCOOH': 46.03,    # муравьиная кислота
    'C2H5COOH': 74.08, # пропионовая кислота
    'CH3COOH': 60.05   # уксусная кислота
}

KA_VALUES = {
    'HCOOH': 10**-3.75,   # муравьиная
    'C2H5COOH': 10**-4.87, # пропионовая
    'CH3COOH': 10**-4.76   # уксусная
}

KW = 1e-14

# === Вспомогательные функции ===
def get_positive_float(prompt: str) -> float:
    """Получает неотрицательное число от пользователя."""
    while True:
        try:
            value = float(input(prompt))
            if value < 0:
                print("Значение не может быть отрицательным!")
                continue
            return value
        except ValueError:
            print("Введите корректное число!")

def calculate_base_conc(c_ha: float, c_a: float, ka: float, h: float) -> float:
    """Вычисляет равновесную концентрацию сопряжённого основания [A⁻]."""
    if c_ha == 0 and c_a == 0:
        return 0.0
    return (ka * c_ha + c_a * h) / (h + ka)

def charge_balance(h: float, c_ha_list, c_a_list, ka_list, c_na: float) -> float:
    """Уравнение электронейтральности: [H⁺] + [Na⁺] - ([OH⁻] + Σ[A⁻]) = 0"""
    oh = KW / max(h, 1e-16)
    total_anions = sum(
        calculate_base_conc(c_ha, c_a, ka, h)
        for c_ha, c_a, ka in zip(c_ha_list, c_a_list, ka_list)
    )
    return h + c_na - (oh + total_anions)

# === Основная логика ===
def main():
    print("=== Расчёт pH смеси слабых кислот ===\n")
    
    # Ввод данных
    mass_hcoona = get_positive_float("Масса формиата натрия (г): ")
    mass_hcooh = get_positive_float("Масса муравьиной кислоты (г): ")
    mass_c2h5cooh = get_positive_float("Масса пропионовой кислоты (г): ")
    mass_ch3cooh = get_positive_float("Масса уксусной кислоты (г): ")
    volume = get_positive_float("Объём раствора (л): ")

    if volume == 0:
        print("Ошибка: объём не может быть нулевым!")
        return

    # Расчёт концентраций
    c_hcoona = mass_hcoona / MOLAR_MASSES['HCOONa'] / volume
    c_hcooh = mass_hcooh / MOLAR_MASSES['HCOOH'] / volume
    c_c2h5cooh = mass_c2h5cooh / MOLAR_MASSES['C2H5COOH'] / volume
    c_ch3cooh = mass_ch3cooh / MOLAR_MASSES['CH3COOH'] / volume

    # Списки для обобщённого расчёта
    c_ha_list = [c_hcooh, c_c2h5cooh, c_ch3cooh]  # кислоты
    c_a_list = [c_hcoona, 0.0, 0.0]                # сопряжённые основания
    ka_list = [KA_VALUES['HCOOH'], KA_VALUES['C2H5COOH'], KA_VALUES['CH3COOH']]

    c_na = c_hcoona  # Na⁺ только от формиата

    # Метод Ньютона
    initial_guesses = [1e-2, 1e-3, 1e-4]
    best_pH = None
    best_balance = float('inf')

    for h0 in initial_guesses:
        h = h0
        for _ in range(100):
            f = charge_balance(h, c_ha_list, c_a_list, ka_list, c_na)
            if abs(f) < 1e-10:
                pH = -math.log10(max(h, 1e-16))
                if abs(f) < abs(best_balance):
                    best_pH = pH
                    best_balance = f
                break
            # Численная производная
            delta = 1e-8
            f_prime = (charge_balance(h + delta, c_ha_list, c_a_list, ka_list, c_na) - f) / delta
            if abs(f_prime) < 1e-12:
                break
            h_new = h - f / f_prime
            if h_new <= 0:
                h_new = 1e-10
            if abs(h_new - h) < 1e-12:
                pH = -math.log10(max(h_new, 1e-16))
                if abs(f) < abs(best_balance):
                    best_pH = pH
                    best_balance = f
                break
            h = h_new

    # Вывод
    if best_pH is not None:
        print(f"\nРезультат: pH = {best_pH:.3f}")
        print(f"Баланс заряда: {best_balance:.2e}")
    else:
        print("\nНе удалось рассчитать pH. Проверьте входные данные.")

if __name__ == "__main__":
    main()
