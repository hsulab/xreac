"""ReaxFF on caller-supplied directed, image-resolved neighbor arrays.

Uses the same equations and units as energy.py, retaining the dense replicated
model as an independent reference. Atom reductions sum all neighbor images;
angles and torsions identify an atom by its index AND accumulated lattice shift.
LAMMPS/PuReMD equation attribution and licensing: see energy.py and NOTICE.
"""
from itertools import combinations

import autograd.numpy as np
from autograd.tracer import getval
import numpy as onp

from .energy import (BOND_CUT, HBOND_CUT, HBOND_THRESHOLD, BOND_GRAPH_CUT,
                     THB_CUT, THB_PRODUCT_CUT, COMPONENTS, C_ELE, QEQ_COULOMB,
                     SELF_CONVERSION, positive_power)
from .neighbors import Neighbors


class NeighborEnergyModel:
    def __init__(self, ff, symbols, neighbors, cell=None, pbc=None):
        ff.validate_model(symbols)
        self.ff, self.symbols, self.n = ff, tuple(symbols), len(symbols)
        self.g, self.vdw_type = ff.general, ff.vdw_type
        self.edges = Neighbors(neighbors, self.n, cell, pbc)
        self.i, self.j = self.edges.i, self.edges.j
        self.a = {k: onp.array([ff.atoms[s][k] for s in symbols]) for k in ff.atoms[symbols[0]]}
        self.p = {k: onp.array([ff.pairs[symbols[i], symbols[j]][k]
                               for i, j in zip(self.i, self.j)])
                  for k in ff.pairs[symbols[0], symbols[0]]}

    def electrostatics(self, r, fixed_charges=None):
        u = np.minimum(r/self.g[12], 1.)
        taper = np.where(r < self.g[12], 1+u**4*(-35+u*(84+u*(-70+20*u))), 0.)
        shield = taper/(r**3+self.p["gamma"])**(1/3)
        if fixed_charges is not None:
            return fixed_charges, taper, shield
        # Nonzero-shift self edges contribute to the QEq diagonal. There is
        # exactly one charge degree of freedom per input atom, even in tiny cells.
        h = QEQ_COULOMB*self.edges.pair_sum(shield)+np.diag(self.a["eta"])
        ones = np.ones((self.n, 1))
        kkt = np.concatenate((np.concatenate((h, ones), axis=1),
                              np.concatenate((ones.T, np.zeros((1, 1))), axis=1)), axis=0)
        q = np.linalg.solve(kkt, np.concatenate((-self.a["chi"], np.zeros(1))))[:self.n]
        return q, taper, shield

    def bond_orders(self, r):
        p, a, g, i, j = self.p, self.a, self.g, self.i, self.j
        cutoff = .01*g[29]
        sigma_ok = (a["r_s"][i] > 0) & (a["r_s"][j] > 0)
        sigma = np.where(sigma_ok, (1+cutoff)*np.exp(p["p_bo1"]*
                         (r/np.where(sigma_ok, p["r_s"], 1))**p["p_bo2"]), 0)
        pi_ok = (a["r_pi"][i] > 0) & (a["r_pi"][j] > 0)
        pp_ok = (a["r_pi_pi"][i] > 0) & (a["r_pi_pi"][j] > 0)
        pi = np.where(pi_ok, np.exp(p["p_bo3"]*(r/np.where(pi_ok, p["r_p"], 1))**p["p_bo4"]), 0)
        pp = np.where(pp_ok, np.exp(p["p_bo5"]*(r/np.where(pp_ok, p["r_pp"], 1))**p["p_bo6"]), 0)
        raw = sigma+pi+pp
        mask = (getval(raw) >= cutoff) & (getval(r) <= min(BOND_CUT, g[12]))
        bo = np.where(mask, raw-cutoff, 0)
        pi, pp = np.where(mask, pi, 0), np.where(mask, pp, 0)
        total = self.edges.atom_sum(bo)
        d, db = total-a["valency"], total-a["valency_boc"]
        f2 = np.exp(-g[0]*d[i])+np.exp(-g[0]*d[j])
        f3 = -np.log(.5*(np.exp(-g[1]*d[i])+np.exp(-g[1]*d[j])))/g[1]
        vi, vj = a["valency"][i], a["valency"][j]
        f1 = np.where(p["ovc"] >= .001, .5*((vi+f2)/(vi+f2+f3)+(vj+f2)/(vj+f2+f3)), 1)
        f4 = 1/(1+np.exp(-(p["p_boc4"]*bo**2-db[i])*p["p_boc3"]+p["p_boc5"]))
        f5 = 1/(1+np.exp(-(p["p_boc4"]*bo**2-db[j])*p["p_boc3"]+p["p_boc5"]))
        corr = f1*np.where(p["v13cor"] >= .001, f4*f5, 1)
        bo, pi, pp = bo*corr, pi*corr*f1, pp*corr*f1
        return tuple(np.where(v >= 1e-10, v, 0.) for v in (bo, bo-pi-pp, pi, pp))

    def properties(self, x, charges):
        _, r = self.edges.geometry(x)
        bo = self.bond_orders(r)[0]
        total = self.edges.atom_sum(bo)
        de = total-self.a["valency_e"]
        half = onp.trunc(getval(de)/2)
        nlp = np.exp(-self.g[15]*(2+de-2*half)**2)-half
        center = np.sum(x*self.a["mass"][:, None], axis=0)/np.sum(self.a["mass"])
        return dict(bond_orders=self.edges.pair_sum(bo), total_bond_orders=total,
                    lone_pairs=nlp, bond_counts=self.edges.atom_sum(
                        (getval(bo) > BOND_GRAPH_CUT).astype(float)).astype(int),
                    dipole=np.sum((x-center)*charges[:, None], axis=0))

    def components(self, x, fixed_charges=None):
        p, a, g, i, j = self.p, self.a, self.g, self.i, self.j
        vectors, r = self.edges.geometry(x)
        q, taper, shield = self.electrostatics(r, fixed_charges)
        bo, sigma, pi, pp = self.bond_orders(r)
        total = self.edges.atom_sum(bo)
        d, db = total-a["valency"], total-a["valency_boc"]
        de = total-a["valency_e"]
        half = onp.trunc(getval(de)/2)
        vlpex = de-2*half
        nlp = np.exp(-g[15]*(2+vlpex)**2)-half
        dlp = .5*(a["valency_e"]-a["valency"])-nlp
        dlpt = np.where(a["mass"] > 21, 0, dlp)
        e = {k: 0. for k in COMPONENTS}
        e["bond"] = -.5*np.sum(p["De_s"]*sigma*np.exp(p["p_be1"]*(1-positive_power(sigma, p["p_be2"])))
                              +p["De_p"]*pi+p["De_pp"]*pp)
        e["lone_pair"] = np.sum(a["p_lp2"]*dlp/(1+np.exp(-75*dlp)))
        self.special_bond_corrections(e, bo, total, d)
        dfvl = (a["mass"] <= 21).astype(float)
        ov1 = self.edges.atom_sum(p["p_ovun1"]*p["De_s"]*bo)
        ov2 = self.edges.atom_sum((d[j]-dfvl[i]*dlpt[j])*(pi+pp))
        dc = d-dfvl*dlpt/(1+g[32]*np.exp(g[31]*ov2))
        over = ov1*dc/(dc+a["valency"]+1e-8)/(1+np.exp(a["p_ovun2"]*dc))
        under = -a["p_ovun5"]*(1-np.exp(g[6]*dc))/(1+np.exp(-a["p_ovun2"]*dc))/(1+g[8]*np.exp(g[9]*ov2))
        e["atom"] = np.sum(over+under)
        f13 = (r**g[28]+(1/p["gamma_w"])**g[28])**(1/g[28]) if self.vdw_type in (1, 3) else r
        ev = np.exp(.5*p["alpha"]*(1-f13/p["r_vdW"]))
        e["vdw"] = .5*np.sum(taper*p["D"]*(ev**2-2*ev))
        if self.vdw_type in (2, 3):
            e["vdw"] += .5*np.sum(taper*p["ecore"]*np.exp(p["acore"]*(1-r/p["rcore"])))
        e["coulomb"] = .5*C_ELE*np.sum(q[i]*q[j]*shield)
        e["qeq"] = SELF_CONVERSION*np.sum(a["chi"]*q+.5*a["eta"]*q*q)
        neighbors = [row[getval(bo)[row] > THB_CUT] for row in self.edges.rows]
        self.angles(e, vectors, r, bo, pi, pp, total, d, db, nlp, vlpex, neighbors)
        self.torsions(e, vectors, r, bo, pi, db, neighbors)
        self.hydrogen_bonds(e, vectors, r, bo)
        return np.stack([e[k] for k in COMPONENTS]), q

    def special_bond_corrections(self, e, bo, total, d):
        a, g, i, j = self.a, self.g, self.i, self.j
        carbon = onp.array([s.upper() == "C" for s in self.symbols])
        if g[5] > .001 and onp.any(carbon):
            vov3 = bo-d[i]-.040*d[i]**4
            mask = carbon[i] & carbon[j] & (getval(vov3) > 3)
            e["lone_pair"] += np.sum(np.where(mask, g[5]*(vov3-3)**2, 0))
        masses = a["mass"]
        co = ((masses[i] == 12.) & (masses[j] == 15.999)) | ((masses[j] == 12.) & (masses[i] == 15.999))
        allowed = co if int(g[37]) != 2 else onp.ones(len(i), dtype=bool)
        selected = onp.flatnonzero(allowed & (getval(bo) >= 1.) & self.edges.half)
        if len(selected):
            i, j, b = i[selected], j[selected], bo[selected]
            stabilization = g[10]*np.exp(-g[7]*(b-2.5)**2)/(1+25*np.exp(g[4]*(d[i]+d[j])))
            stabilization *= np.exp(-g[3]*(total[i]-b))+np.exp(-g[3]*(total[j]-b))
            e["bond"] += np.sum(stabilization)

    def hydrogen_bonds(self, e, vectors, r, bo):
        roles, bv, rv = self.a["p_hbond"], getval(bo), getval(r)
        pairs, params = [], []
        for j in onp.flatnonzero(roles == 1):
            row = self.edges.rows[j]
            heavy = row[roles[self.j[row]] == 2]
            donors = heavy[bv[heavy] >= HBOND_THRESHOLD]
            acceptors = heavy[rv[heavy] <= min(HBOND_CUT, self.g[12])]
            for donor in donors:
                for acceptor in acceptors:
                    # Exclude only the SAME image; other images of the donor
                    # remain valid acceptors even when their atom indices match.
                    if donor == acceptor:
                        continue
                    key = (self.symbols[self.j[donor]], self.symbols[j], self.symbols[self.j[acceptor]])
                    prm = self.ff.hydrogen_bonds.get(key)
                    if prm is not None and prm[0] > 0:
                        pairs.append((donor, acceptor))
                        params.append(prm)
        if not pairs:
            return
        donor, acceptor = onp.array(pairs).T
        r0, hb1, hb2, hb3 = onp.array(params).T
        distance = r[acceptor]
        cos = np.sum(vectors[donor]*vectors[acceptor], axis=1)/(r[donor]*distance)
        sin4 = .25*(1-np.clip(cos, -1., 1.))**2
        e["hydrogen_bond"] = np.sum(hb1*(1-np.exp(-hb2*bo[donor]))
            *np.exp(-hb3*(r0/distance+distance/r0-2))*sin4)

    def angles(self, e, vectors, r, bo, pi, pp, total, d, db, nlp, vlpex, neighbors):
        a, g, bv = self.a, self.g, getval(bo)
        pairs, params = [], []
        for j, row in enumerate(neighbors):
            for left, right in combinations(row, 2):
                if bv[left]*bv[right] <= THB_PRODUCT_CUT:
                    continue
                key = (self.symbols[self.j[left]], self.symbols[j], self.symbols[self.j[right]])
                for prm in self.ff.angles.get(key, []):
                    if abs(prm[1]) > .001:
                        pairs.append((left, right))
                        params.append(prm)
        if not pairs:
            return
        left, right = onp.array(pairs).T
        i, j, k = self.j[left], self.i[left], self.j[right]
        theta00, v1, v2, coa1, v7, pen1, v4 = onp.array(params).T
        prod = np.exp(-self.edges.atom_sum(bo**8))
        sbo = self.edges.atom_sum(pi+pp)+(1-prod)*(-db-g[33]*np.where(getval(vlpex) >= 0, 0, nlp))
        sbo2 = np.where(sbo <= 0, 0, np.where(sbo <= 1, positive_power(sbo, g[16]),
                        np.where(sbo < 2, 2-positive_power(2-sbo, g[16]), 2)))
        theta0 = (180-theta00*(1-np.exp(-g[17]*(2-sbo2[j]))))*(3.14159265/180)
        cos = np.sum(vectors[left]*vectors[right], axis=1)/(r[left]*r[right])
        theta = np.arccos(np.clip(cos, -1+1e-14, 1-1e-14))
        bij, bjk = bo[left]-THB_CUT, bo[right]-THB_CUT
        f7 = (1-np.exp(-a["p_val3"][j]*bij**v4))*(1-np.exp(-a["p_val3"][j]*bjk**v4))
        ex6, ex7 = np.exp(g[14]*db[j]), np.exp(-v7*db[j])
        f8 = a["p_val5"][j]-(a["p_val5"][j]-1)*(2+ex6)/(1+ex6+ex7)
        exangle = np.exp(-v2*(theta0-theta)**2)
        e["angle"] = np.sum(f7*f8*v1*np.where(v1 >= 0, 1-exangle, -exangle))
        ex3, ex4 = np.exp(-g[20]*d[j]), np.exp(g[21]*d[j])
        e["penalty"] = np.sum(pen1*(2+ex3)/(1+ex3+ex4)*np.exp(-g[19]*((bij-2)**2+(bjk-2)**2)))
        e["angle_conjugation"] = np.sum(coa1/(1+np.exp(g[2]*(total[j]-a["valency_val"][j])))
            *np.exp(-g[38]*((total[i]-bij)**2+(total[k]-bjk)**2))
            *np.exp(-g[30]*((bij-1.5)**2+(bjk-1.5)**2)))

    def torsions(self, e, vectors, r, bo, pi, db, neighbors):
        chains, params, bv = [], [], getval(bo)
        shifts = self.edges.shifts
        for middle in onp.flatnonzero(self.edges.half & (bv > THB_CUT)):
            j, k = self.i[middle], self.j[middle]
            for left in neighbors[j]:
                if left == middle:
                    continue
                for right in neighbors[k]:
                    if right == self.edges.reverse[middle] or bv[left]*bv[middle]*bv[right] <= THB_CUT:
                        continue
                    i, l = self.j[left], self.j[right]
                    if i == l and onp.array_equal(shifts[left], shifts[middle]+shifts[right]):
                        continue
                    prm = self.ff.torsions.get(tuple(self.symbols[t] for t in (i, j, k, l)))
                    if prm is not None and onp.any(prm[[0, 1, 2, 4]] != 0):
                        chains.append((left, middle, right))
                        params.append(prm)
        if not chains:
            return
        left, middle, right = onp.array(chains).T
        j, k = self.i[middle], self.j[middle]
        v1, v2, v3, tor1, cot1 = onp.array(params).T
        u, v, w = vectors[left], vectors[middle], vectors[right]
        normal1, normal2 = np.cross(u, v), np.cross(w, v)
        area1, area2 = np.sum(normal1**2, axis=1), np.sum(normal2**2, axis=1)
        sine1_sq = area1/(r[left]*r[middle])**2
        sine2_sq = area2/(r[middle]*r[right])**2
        if onp.any(getval(sine1_sq) < 1e-20) or onp.any(getval(sine2_sq) < 1e-20):
            raise ValueError("Collinear atoms in an active torsion: the dihedral derivative is undefined; perturb the geometry")
        sinprod = np.sqrt(sine1_sq*sine2_sq)
        cw = np.sum(normal1*normal2, axis=1)/np.sqrt(area1*area2)
        b1, b2, b3 = bo[left]-THB_CUT, bo[middle]-THB_CUT, bo[right]-THB_CUT
        ds = db[j]+db[k]
        ex3, ex4 = np.exp(-self.g[24]*ds), np.exp(self.g[25]*ds)
        f11 = (2+ex3)/(1+ex3+ex4)
        ex1 = np.exp(tor1*(2-pi[middle]-f11)**2)
        fn10 = (1-np.exp(-self.g[23]*b1))*(1-np.exp(-self.g[23]*b2))*(1-np.exp(-self.g[23]*b3))
        cv = .5*(v1*(1+cw)+v2*ex1*(2-2*cw*cw)+v3*(1+4*cw**3-3*cw))
        e["torsion"] = np.sum(fn10*sinprod*cv)
        fn12 = np.exp(-self.g[27]*((b1-1.5)**2+(b2-1.5)**2+(b3-1.5)**2))
        e["conjugation"] = np.sum(cot1*fn12*(1+(cw*cw-1)*sinprod))
