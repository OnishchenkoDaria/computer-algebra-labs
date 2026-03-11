from itertools import product
from math import isqrt
from sympy import symbols, Poly, ZZ, div, interpolate, expand

x = symbols('x')


def integer_divisors(n: int):
    """Повертає всі цілі дільники числа n."""
    if n == 0:
        raise ValueError("Для 0 множина дільників нескінченна.")
    n_abs = abs(n)
    divs = set()
    for d in range(1, isqrt(n_abs) + 1):
        if n_abs % d == 0:
            divs.add(d)
            divs.add(-d)
            divs.add(n_abs // d)
            divs.add(-(n_abs // d))
    return sorted(divs)


def poly_from_values(points, values):
    """
    Відновлює поліном за значеннями у точках через інтерполяцію.
    points: список x_i
    values: список y_i
    """
    expr = interpolate(list(zip(points, values)), x)
    return Poly(expand(expr), x, domain='QQ')


def has_integer_coeffs(poly: Poly) -> bool:
    """Перевіряє, що всі коефіцієнти полінома цілі."""
    return all(c.is_integer for c in poly.all_coeffs())


def primitive_integer_poly(poly: Poly) -> Poly:
    """
    Переводить поліном з QQ у Z, якщо це можливо через очищення знаменників.
    Для методу Кронекера краще залишати лише цілочисельні кандидати,
    тому тут тільки контрольна нормалізація.
    """
    coeffs = poly.all_coeffs()
    if all(c.is_integer for c in coeffs):
        return Poly(poly.as_expr(), x, domain=ZZ)
    return poly


def kronecker_find_factor(f: Poly):
    """
    Знаходить один нетривіальний дільник полінома f методом Кронекера
    або повертає None, якщо не знайшов.
    """
    f = Poly(f, x, domain=ZZ)
    n = f.degree()

    if n <= 1:
        return None

    # 1. Шукаємо цілі корені серед 0, 1, ..., floor(n/2) та їх заперечень
    # Для класичного опису часто перевіряють f(i)=0 на невеликому наборі точок,
    # але на практиці корисно перевірити дільники вільного члена.
    c0 = f.eval(0)
    if c0 != 0:
        for r in integer_divisors(int(c0)):
            if f.eval(r) == 0:
                return Poly(x - r, x, domain=ZZ)

    # 2. Шукаємо дільник степеня m = 1..floor(n/2)
    max_deg = n // 2

    for m in range(1, max_deg + 1):
        # Беремо m+1 точок
        points = list(range(m + 1))

        # Якщо в якійсь точці f(i)=0, краще змістити набір точок
        shift = 0
        while True:
            shifted_points = [p + shift for p in points]
            vals = [f.eval(p) for p in shifted_points]
            if all(v != 0 for v in vals):
                break
            shift += 1

        divisor_lists = [integer_divisors(int(v)) for v in vals]

        for candidate_values in product(*divisor_lists):
            g_q = poly_from_values(shifted_points, candidate_values)

            # Степінь не має перевищувати m і не має бути 0
            if g_q.degree() <= 0 or g_q.degree() > m:
                continue

            # Залишаємо лише цілочисельні кандидати
            if not has_integer_coeffs(g_q):
                continue

            g = primitive_integer_poly(g_q)

            # Відкидаємо ±1
            if g.degree() == 0:
                continue

            # Перевіряємо точне ділення
            q, r = div(f, g, domain=ZZ)
            if r.is_zero:
                # Не беремо тривіальний випадок g = ±f
                if g.degree() < f.degree():
                    return g

    return None


def kronecker_factorization(f: Poly):
    """
    Повний рекурсивний розклад методом Кронекера.
    Повертає список нерозкладних множників над Z[x].
    """
    f = Poly(f, x, domain=ZZ)

    if f.degree() <= 1:
        return [f]

    factor = kronecker_find_factor(f)
    if factor is None:
        return [f]

    q, r = div(f, factor, domain=ZZ)
    if not r.is_zero:
        return [f]

    return kronecker_factorization(factor) + kronecker_factorization(q)


if __name__ == "__main__":
    f = Poly(x**5 - x**4 - 2*x**3 - 8*x**2 + 6*x - 1, x, domain=ZZ)

    factors = kronecker_factorization(f)

    print("Поліном:")
    print(f.as_expr())
    print("\nЗнайдені множники:")
    for fac in factors:
        print(fac.as_expr())