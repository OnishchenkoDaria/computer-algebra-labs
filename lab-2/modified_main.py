from sympy import symbols, Poly, QQ, expand, groebner
from collections import Counter

x, y, z = symbols('x y z')
GENS = (x, y, z)   # lex: x > y > z

def P(expr):
    return Poly(expand(expr), *GENS, domain=QQ)

# Початкова система
f1 = P(x*y - x*z + y**2)
f2 = P(y*z - x**2 + x**2*y)
f3 = P(x - x*y + y)

def monomial_to_expr(m):
    result = 1
    for var, power in zip(GENS, m):
        result *= var**power
    return result

def monomial_divides(a, b):
    return all(ai <= bi for ai, bi in zip(a, b))

def monomial_quotient(num, den):
    return tuple(ni - di for ni, di in zip(num, den))

def monomial_lcm(a, b):
    return tuple(max(ai, bi) for ai, bi in zip(a, b))

def relatively_prime(m1, m2):
    # product criterion
    return all(min(a, b) == 0 for a, b in zip(m1, m2))

def LM(p):
    mons = p.monoms()
    return None if not mons else mons[0]

def LC(p):
    coeffs = p.coeffs()
    return 0 if not coeffs else coeffs[0]

def LT_expr(p):
    if p.as_expr() == 0:
        return 0
    return LC(p) * monomial_to_expr(LM(p))

def s_polynomial(p, q):
    mp = LM(p)
    mq = LM(q)
    l = monomial_lcm(mp, mq)

    mult_p = monomial_to_expr(monomial_quotient(l, mp))
    mult_q = monomial_to_expr(monomial_quotient(l, mq))

    sp = mult_p / LC(p) * p.as_expr() - mult_q / LC(q) * q.as_expr()
    return P(sp)

def reduce_poly(f, G, trace=False):
    p = P(f.as_expr() if isinstance(f, Poly) else f)
    r = P(0)
    steps = []

    while p.as_expr() != 0:
        reduced = False
        mp = LM(p)
        cp = LC(p)

        for g in G:
            if g.as_expr() == 0:
                continue

            mg = LM(g)
            cg = LC(g)

            if monomial_divides(mg, mp):
                qmon = monomial_quotient(mp, mg)
                factor = cp / cg * monomial_to_expr(qmon)
                new_p = P(p.as_expr() - factor * g.as_expr())

                if trace:
                    steps.append({
                        "type": "reduce",
                        "current": p.as_expr(),
                        "by": g.as_expr(),
                        "factor": factor,
                        "result": new_p.as_expr()
                    })

                p = new_p
                reduced = True
                break

        if not reduced:
            lt = LT_expr(p)
            r = P(r.as_expr() + lt)
            new_p = P(p.as_expr() - lt)

            if trace:
                steps.append({
                    "type": "move_to_remainder",
                    "current": p.as_expr(),
                    "lt": lt,
                    "new_current": new_p.as_expr(),
                    "remainder": r.as_expr()
                })

            p = new_p

    return r, steps

def make_monic(p):
    if p.as_expr() == 0:
        return p
    return P(p.as_expr() / LC(p))

def minimal_groebner_basis(G):
    Gm = [make_monic(g) for g in G if g.as_expr() != 0]

    result = []
    for i, gi in enumerate(Gm):
        keep = True
        for j, gj in enumerate(Gm):
            if i != j and monomial_divides(LM(gj), LM(gi)):
                keep = False
                break
        if keep:
            result.append(gi)

    unique = []
    seen = set()
    for g in result:
        expr = expand(g.as_expr())
        if expr not in seen:
            seen.add(expr)
            unique.append(P(expr))
    return unique

def reduced_groebner_basis(G):
    Gmin = minimal_groebner_basis(G)
    Gred = []

    for i, g in enumerate(Gmin):
        others = [Gmin[j] for j in range(len(Gmin)) if j != i]
        r, _ = reduce_poly(g, others, trace=False)
        Gred.append(make_monic(r))

    final = []
    seen = set()
    for g in Gred:
        if g.as_expr() != 0:
            expr = expand(g.as_expr())
            if expr not in seen:
                seen.add(expr)
                final.append(P(expr))
    return final

# ---------------- БАЗОВИЙ ----------------
def buchberger_basic(F):
    G = [P(f.as_expr() if isinstance(f, Poly) else f) for f in F]
    pairs = [(i, j) for i in range(len(G)) for j in range(i)]
    history = []

    while pairs:
        i, j = pairs.pop(0)

        S = s_polynomial(G[i], G[j])
        R, steps = reduce_poly(S, G, trace=True)

        history.append({
            "pair": (i, j),
            "status": "processed",
            "S": S.as_expr(),
            "R": R.as_expr(),
            "added": (R.as_expr() != 0),
            "steps": steps
        })

        if R.as_expr() != 0:
            G.append(R)
            new_idx = len(G) - 1
            for k in range(new_idx):
                pairs.append((new_idx, k))

    return G, history

# ---------------- ОПТИМІЗОВАНИЙ ----------------
def buchberger_optimized(F, use_product_criterion=True, use_signature_filter=True):
    G = [P(f.as_expr() if isinstance(f, Poly) else f) for f in F]
    pairs = [(i, j) for i in range(len(G)) for j in range(i)]
    history = []

    processed_pairs = set()
    known_polys = set(expand(make_monic(g).as_expr()) for g in G)
    s_computations = 0

    while pairs:
        i, j = pairs.pop(0)

        key = tuple(sorted((i, j)))
        if key in processed_pairs:
            continue
        processed_pairs.add(key)

        # Product criterion
        if use_product_criterion and relatively_prime(LM(G[i]), LM(G[j])):
            history.append({
                "pair": (i, j),
                "status": "skip_product_criterion",
                "added": False
            })
            continue

        S = s_polynomial(G[i], G[j])
        s_computations += 1
        R, steps = reduce_poly(S, G, trace=True)
        R = make_monic(R)

        if R.as_expr() == 0:
            history.append({
                "pair": (i, j),
                "status": "reduced_to_zero",
                "S": S.as_expr(),
                "R": 0,
                "added": False,
                "steps": steps
            })
            continue

        exprR = expand(R.as_expr())

        # simplified F5-like duplicate filter
        if use_signature_filter and exprR in known_polys:
            history.append({
                "pair": (i, j),
                "status": "skip_duplicate_polynomial",
                "S": S.as_expr(),
                "R": exprR,
                "added": False,
                "steps": steps
            })
            continue

        known_polys.add(exprR)
        history.append({
            "pair": (i, j),
            "status": "added",
            "S": S.as_expr(),
            "R": exprR,
            "added": True,
            "steps": steps
        })

        G.append(R)
        new_idx = len(G) - 1
        for k in range(new_idx):
            new_key = tuple(sorted((new_idx, k)))
            if new_key not in processed_pairs:
                pairs.append((new_idx, k))

    return G, history, s_computations


if __name__ == "__main__":
    F = [f1, f2, f3]

    print("Початкова система:")
    for i, f in enumerate(F, 1):
        print(f"f{i} =", f.as_expr())

    G_basic, hist_basic = buchberger_basic(F)
    G_basic_red = reduced_groebner_basis(G_basic)

    basic_s = len(hist_basic)
    basic_added = sum(1 for h in hist_basic if h["added"])

    print("\n=== БАЗОВИЙ АЛГОРИТМ ===")
    print("Кількість пар / обчислень S-поліномів:", basic_s)
    print("Кількість доданих поліномів:", basic_added)
    for i, g in enumerate(G_basic_red, 1):
        print(f"r{i} =", g.as_expr())

    G_opt, hist_opt, opt_s = buchberger_optimized(
        F,
        use_product_criterion=True,
        use_signature_filter=True
    )
    G_opt_red = reduced_groebner_basis(G_opt)

    opt_added = sum(1 for h in hist_opt if h["added"])
    stat = Counter(h["status"] for h in hist_opt)

    print("\n=== ОПТИМІЗОВАНИЙ АЛГОРИТМ ===")
    print("Кількість реально обчислених S-поліномів:", opt_s)
    print("Кількість доданих поліномів:", opt_added)
    print("Статистика:", dict(stat))
    for i, g in enumerate(G_opt_red, 1):
        print(f"r{i} =", g.as_expr())

    print("\n=== ПОРІВНЯННЯ ===")
    print("Базовий:", basic_s)
    print("Оптимізований:", opt_s)
    if basic_s:
        print("Зменшення: {:.2f}%".format((basic_s - opt_s) / basic_s * 100))

    print("\n=== ПЕРЕВІРКА БІБЛІОТЕКОЮ ===")
    print(groebner([f1.as_expr(), f2.as_expr(), f3.as_expr()], x, y, z, order="lex"))