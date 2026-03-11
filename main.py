from itertools import product
from math import isqrt
from sympy import symbols, Poly, ZZ, QQ, div, interpolate, expand

x = symbols('x')


def divisors(n):
    """Пошук усіх цілих дільників числа"""
    n = abs(n)
    result = set()

    for i in range(1, isqrt(n) + 1):
        if n % i == 0:
            result.update({i, -i, n // i, -(n // i)})

    return sorted(result)


def interpolate_poly(points, values):
    """Побудова полінома за значеннями (інтерполяція)"""
    expr = interpolate(list(zip(points, values)), x)
    return Poly(expand(expr), x, domain=QQ)


def kronecker_factor(f):

    f = Poly(f, x, domain=ZZ)

    print("\n--- Початок алгоритму Кронекера ---")
    print("Поліном f(x):", f.as_expr())

    n = f.degree()
    print("Степінь полінома:", n)

    if n <= 1:
        return None

    # межа степеня можливого множника
    max_deg = n // 2
    print("Максимальний степінь шуканого множника:", max_deg)

    # перебір можливих степенів множника
    for m in range(1, max_deg + 1):

        print("\nЕтап: перевірка множників степеня", m)

        points = list(range(m + 1))
        print("Вибрані точки:", points)

        values = [f.eval(i) for i in points]
        print("Значення f(i):", values)

        # пошук дільників
        div_lists = []
        for p, v in zip(points, values):
            d = divisors(int(v))
            div_lists.append(d)
            print(f"Дільники числа f({p}) = {v}:", d)

        # перебір можливих значень g(i)
        for combo in product(*div_lists):

            print("\nПеревірка набору значень g(i):", combo)

            # побудова кандидата g(x)
            g = interpolate_poly(points, combo)

            print("Кандидатний поліном g(x):", g.as_expr())

            # перевірка степеня
            if g.degree() <= 0 or g.degree() > m:
                print("Відкинуто: некоректний степінь")
                continue

            # перевірка цілочисельних коефіцієнтів
            if not all(c.is_integer for c in g.all_coeffs()):
                print("Відкинуто: коефіцієнти не цілі")
                continue

            g = Poly(g.as_expr(), x, domain=ZZ)

            # перевірка ділення
            q, r = div(f, g)

            print("Остача від ділення:", r.as_expr())

            if r.is_zero and g.degree() < f.degree():
                print("Знайдено множник:", g.as_expr())
                return g

    print("Множник не знайдено")
    return None


def kronecker_factorization(f):

    f = Poly(f, x, domain=ZZ)

    factor = kronecker_factor(f)

    if factor is None:
        print("\nПоліном нерозкладний:", f.as_expr())
        return [f]

    q, _ = div(f, factor)

    print("\nРозклад:")
    print("f(x) =", factor.as_expr(), "*", q.as_expr())

    return (
        kronecker_factorization(factor)
        + kronecker_factorization(q)
    )


# тестовий приклад
if __name__ == "__main__":

    print("Номер студента у списку - 18 : варіант 4\n")
    print("\n==================== TASK 18 ====================\n")
    f = Poly(x**5 + x**4 - 21*x**3 - 45*x**2, x)

    factors = kronecker_factorization(f)

    print("\n--- Підсумкові множники ---")
    for g in factors:
        print(g.as_expr())

    print("\n==================== TASK 24 ====================\n")
    f = Poly(2 * x ** 4 + 14 * x ** 3 + 12 * x ** 2 - 56 * x - 80, x)

    factors = kronecker_factorization(f)

    print("\n--- Підсумкові множники ---")
    for g in factors:
        print(g.as_expr())