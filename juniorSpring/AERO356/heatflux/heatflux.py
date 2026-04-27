import math

SIGMA = 5.67e-8  # W/m^2/K^4


class SpacecraftThermal:
    def __init__(self, area, emissivity, name="Spacecraft"):
        if area <= 0:
            raise ValueError("Area must be positive.")
        if emissivity <= 0:
            raise ValueError("Emissivity must be positive.")

        self.name = name
        self.A = float(area)          # radiating area [m^2]
        self.eps = float(emissivity)  # emissivity [-]

        self.inputs = {}              # heat inputs [W]
        self.outputs = {}             # non-radiative outputs [W]

    # =========================================
    # INTERNAL HELPERS
    # =========================================

    def _add_term(self, storage, key, Q, category, notes=""):
        if Q < 0:
            raise ValueError("Flux must be nonnegative.")
        storage[key] = {
            "Q": float(Q),
            "category": category,
            "notes": notes
        }

    def _sum_terms(self, storage):
        return sum(term["Q"] for term in storage.values())

    # =========================================
    # INPUT FLUX TERMS
    # =========================================

    def addSolar(self, *args, key="solar", notes=""):
        """
        Options:
            sc.addSolar(Q)
            sc.addSolar(alpha, S, A_proj)

        Computes:
            Q_solar = alpha * S * A_proj
        """
        if len(args) == 1:
            Q = args[0]
        elif len(args) == 3:
            alpha, S, A_proj = args
            Q = alpha * S * A_proj
        else:
            raise ValueError("addSolar expects either (Q) or (alpha, S, A_proj)")

        self._add_term(self.inputs, key, Q, "solar", notes)

    def addAlbedo(self, *args, key="albedo", notes=""):
        """
        Options:
            sc.addAlbedo(Q)
            sc.addAlbedo(alpha, S, Ab, A_proj)
            sc.addAlbedo(alpha, S, Ab, A_proj, G, theta)

        theta must be in radians.

        Computes:
            Q_albedo = alpha * S * Ab * A_proj * G * cos(theta)
        """
        if len(args) == 1:
            Q = args[0]
        elif len(args) == 4:
            alpha, S, Ab, A_proj = args
            G = 1.0
            theta = 0.0
            Q = alpha * S * Ab * A_proj * G * math.cos(theta)
        elif len(args) == 6:
            alpha, S, Ab, A_proj, G, theta = args
            Q = alpha * S * Ab * A_proj * G * math.cos(theta)
        else:
            raise ValueError(
                "addAlbedo expects either (Q), (alpha, S, Ab, A_proj), "
                "or (alpha, S, Ab, A_proj, G, theta)"
            )

        self._add_term(self.inputs, key, Q, "albedo", notes)

    def addEarthIR(self, *args, key="earth_ir", notes=""):
        """
        Options:
            sc.addEarthIR(Q)
            sc.addEarthIR(eps_sc, eps_earth, T_earth, view_factor)
            sc.addEarthIR(eps_sc, eps_earth, T_earth, view_factor, area)

        Computes:
            Q_IR = eps_sc * eps_earth * sigma * T_earth^4 * area * F
        """
        if len(args) == 1:
            Q = args[0]
        elif len(args) == 4:
            eps_sc, eps_earth, T_earth, view_factor = args
            Q = eps_sc * eps_earth * SIGMA * T_earth**4 * self.A * view_factor
        elif len(args) == 5:
            eps_sc, eps_earth, T_earth, view_factor, area = args
            Q = eps_sc * eps_earth * SIGMA * T_earth**4 * area * view_factor
        else:
            raise ValueError(
                "addEarthIR expects either (Q), "
                "(eps_sc, eps_earth, T_earth, view_factor), "
                "or (eps_sc, eps_earth, T_earth, view_factor, area)"
            )

        self._add_term(self.inputs, key, Q, "earth_ir", notes)

    def addInternal(self, Q, key="internal", notes=""):
        """
        Internal heat source [W].
        """
        self._add_term(self.inputs, key, Q, "internal", notes)

    def addInput(self, name, Q, notes="", category="custom_input"):
        self._add_term(self.inputs, name, Q, category, notes)

    # =========================================
    # OUTPUT TERMS
    # =========================================

    def addOutput(self, name, Q, notes="", category="custom_output"):
        self._add_term(self.outputs, name, Q, category, notes)

    # =========================================
    # REMOVE / EDIT
    # =========================================

    def removeInput(self, name):
        if name in self.inputs:
            del self.inputs[name]

    def removeOutput(self, name):
        if name in self.outputs:
            del self.outputs[name]

    def updateNotes(self, name, notes, storage="input"):
        if storage == "input":
            if name not in self.inputs:
                raise KeyError(f"No input term named '{name}'")
            self.inputs[name]["notes"] = notes
        elif storage == "output":
            if name not in self.outputs:
                raise KeyError(f"No output term named '{name}'")
            self.outputs[name]["notes"] = notes
        else:
            raise ValueError("storage must be 'input' or 'output'")

    # =========================================
    # THERMAL CALCULATIONS
    # =========================================

    def totalInput(self):
        return self._sum_terms(self.inputs)

    def totalOutput(self):
        return self._sum_terms(self.outputs)

    def totalByCategory(self, category, storage="input"):
        if storage == "input":
            source = self.inputs
        elif storage == "output":
            source = self.outputs
        else:
            raise ValueError("storage must be 'input' or 'output'")

        return sum(term["Q"] for term in source.values() if term["category"] == category)

    def nonRadiativeNetFlux(self):
        """
        Returns:
            Q_in - Q_out

        This does NOT include radiation.
        """
        return self.totalInput() - self.totalOutput()

    def radiativeLoss(self, T):
        """
        Q_rad = eps * sigma * A * T^4
        """
        if T < 0:
            raise ValueError("Temperature must be nonnegative.")
        return self.eps * SIGMA * self.A * T**4

    def netFlux(self, T=None):
        """
        If T is None:
            returns Q_in - Q_out

        If T is provided:
            returns Q_in - Q_out - Q_rad(T)
        """
        q_net = self.nonRadiativeNetFlux()
        if T is not None:
            q_net -= self.radiativeLoss(T)
        return q_net

    def temperatureFromRadiatedFlux(self, Q_rad):
        """
        Solve:
            Q_rad = eps * sigma * A * T^4
        """
        if Q_rad <= 0:
            raise ValueError("Radiated flux must be positive.")
        return (Q_rad / (self.eps * SIGMA * self.A)) ** 0.25

    def equilibriumTemperature(self):
        """
        Standard steady-state equilibrium:
            Q_in - Q_out - Q_rad = 0

        So:
            Q_rad = Q_in - Q_out
        """
        Q_rad = self.nonRadiativeNetFlux()
        return self.temperatureFromRadiatedFlux(Q_rad)

    def temperatureForNetGain(self, Q_gain):
        """
        Solve problems of the form:
            Q_in - Q_out - Q_rad = Q_gain

        Rearranged:
            Q_rad = Q_in - Q_out - Q_gain

        This matches homework statements like:
            'If the heat gain is 5 W, what is the temperature?'

        Example:
            Qsolar + Qalb + QIR - Qrad = 5
            -> Qrad = Qsolar + Qalb + QIR - 5
        """
        if Q_gain < 0:
            raise ValueError("Q_gain must be nonnegative.")

        Q_rad = self.nonRadiativeNetFlux() - Q_gain
        return self.temperatureFromRadiatedFlux(Q_rad)

    def temperatureFromBalance(self, rhs):
        """
        General solver for:
            Q_in - Q_out - Q_rad = rhs

        Then:
            Q_rad = (Q_in - Q_out) - rhs

        Examples:
            rhs = 0   -> standard equilibrium
            rhs = 5   -> net heat gain of 5 W
            rhs = -2  -> net cooling of 2 W
        """
        Q_rad = self.nonRadiativeNetFlux() - rhs
        return self.temperatureFromRadiatedFlux(Q_rad)

    def temperatureAtCurrentState(self):
        return self.equilibriumTemperature()

    def Tmax(self, scale_solar=1.0, rhs=0.0):
        """
        Scales solar input terms and solves:
            Q_in - Q_out - Q_rad = rhs
        """
        if scale_solar < 0:
            raise ValueError("scale_solar must be nonnegative.")

        original_values = {}

        for key, term in self.inputs.items():
            if term["category"] == "solar":
                original_values[key] = term["Q"]
                term["Q"] *= scale_solar

        try:
            return self.temperatureFromBalance(rhs)
        finally:
            for key, value in original_values.items():
                self.inputs[key]["Q"] = value

    # =========================================
    # SUMMARY
    # =========================================

    def summary(self, return_string=False):
        lines = []
        lines.append("=" * 60)
        lines.append(f"Thermal Summary: {self.name}")
        lines.append("=" * 60)
        lines.append(f"Radiating area      : {self.A:.6f} m^2")
        lines.append(f"Emissivity          : {self.eps:.6f}")
        lines.append("")

        lines.append("INPUT TERMS")
        lines.append("-" * 60)
        if self.inputs:
            for name, term in self.inputs.items():
                lines.append(f"{name}")
                lines.append(f"  category : {term['category']}")
                lines.append(f"  flux     : {term['Q']:.6f} W")
                lines.append(f"  notes    : {term['notes'] if term['notes'] else '(none)'}")
        else:
            lines.append("(none)")

        lines.append("")
        lines.append("OUTPUT TERMS")
        lines.append("-" * 60)
        if self.outputs:
            for name, term in self.outputs.items():
                lines.append(f"{name}")
                lines.append(f"  category : {term['category']}")
                lines.append(f"  flux     : {term['Q']:.6f} W")
                lines.append(f"  notes    : {term['notes'] if term['notes'] else '(none)'}")
        else:
            lines.append("(none)")

        lines.append("")
        lines.append("TOTALS")
        lines.append("-" * 60)
        lines.append(f"Total input         : {self.totalInput():.6f} W")
        lines.append(f"Total output        : {self.totalOutput():.6f} W")
        lines.append(f"Net before radiation: {self.nonRadiativeNetFlux():.6f} W")

        try:
            T_eq = self.equilibriumTemperature()
            lines.append(f"Equilibrium T       : {T_eq:.6f} K")
        except ValueError as e:
            lines.append(f"Equilibrium T       : not available ({e})")

        lines.append("=" * 60)

        text = "\n".join(lines)

        if return_string:
            return text
        print(text)


if __name__ == "__main__":
    r = 0.5
    A_proj = math.pi * r**2
    A_surf = 4 * math.pi * r**2

    Se = 1366
    bondAlbedo = 0.3
    tempEarth = 255
    RE = 6378e3
    h = 400e3

    sc = SpacecraftThermal(
        area=A_surf,
        emissivity=0.1,
        name="Spherical Spacecraft"
    )

    sc.addSolar(0.1, Se, A_proj, notes="Direct solar")
    rho = (RE + h) / RE
    G = (2 / 3) * ((2 * rho + rho**-2) - (2 + rho**-2) * math.sqrt(rho**2 - 1))
    sc.addAlbedo(0.1, Se, bondAlbedo, A_proj, G, 0.0, notes="Max albedo")

    F = 0.5 * (1 - math.sqrt(1 - (RE / (RE + h))**2))
    sc.addEarthIR(0.1, 0.9, tempEarth, F, notes="Earth IR")

    sc.summary()

    print("\nStandard equilibrium temperature:")
    print(sc.equilibriumTemperature(), "K")

    print("\nTemperature for required net heat gain of 5 W:")
    print(sc.temperatureForNetGain(5.0), "K")