import pyscript
from pyscript.web import page

import raschii
import numpy as np


@pyscript.when("click", "#generate_wave")
def plot_wave():
    """
    This is called when the user clicks the "Generate Wave" button.

    - Read the user input
    - Check the breaking criteria
    - Generate the wave
    - Show the wave profile using SVG
    - Show the wave properties in output text area
    - Show warnings and errors in the info area

    """
    global current_raschii_wave
    current_raschii_wave = None

    page.find("#raschii p.info").textContent = "Generating wave..."
    page.find("#raschii p.warning").textContent = ""
    page.find("#raschii p.error").textContent = ""
    page.find("#raschii_out").textContent = ""

    # Read the user input
    wave_input = WaveInput()
    wave_input.load_from_page()
    if not wave_input.is_ok:
        page.find("#raschii p.error").textContent = wave_input.error_message
        return

    # Check the breaking criteria
    breaking_errors, breaking_warnings = raschii.check_breaking_criteria(
        height=wave_input.height, depth=wave_input.depth, length=wave_input.length
    )
    if breaking_errors:
        page.find("#raschii p.info").textContent = "Breaking criterion triggered!"
        page.find("#raschii p.error").textContent = breaking_errors
        return

    # Load the wave model class and generate the wave
    WaveModelClass = raschii.get_wave_model(wave_input.wave_model_name)[0]
    wave = WaveModelClass(
        height=wave_input.height,
        depth=wave_input.depth,
        length=wave_input.length,
        N=wave_input.order,
    )

    # Show wave model output
    if not breaking_errors:
        page.find("#raschii p.info").textContent = "Wave generated successfully!"
    page.find("#raschii p.warning").textContent = f"{breaking_warnings}\n\n{wave.warnings}"

    # The wave elevation
    t = 0.0  # Time at which to evaluate the wave
    x = np.linspace(0.0, wave_input.length / 2, 200)
    eta = wave.surface_elevation(x, t)
    Uc = wave.velocity(x[0], eta[0], all_points_wet=True).flatten()
    Ut = wave.velocity(x[-1], eta[-1], all_points_wet=True).flatten()

    outputs = [
        "Summary of results:",
        f"  Surface elevation maximum = {eta.max():.3f}",
        f"  Surface elevation minimum = {eta.min():.3f}",
        f"  Horizontal crest particle speed  = {Uc[0]:.3f}",
        f"  Horizontal trough particle speed = {Ut[0]:.3f}",
        f"  Wave number    = {wave.k:.5f}",
        f"  Wave frequency = {wave.omega:.5f}",
        f"  Wave period    = {wave.T:.2f}",
        f"  Phase speed    = {wave.c:.3f}",
        "",
        "Numerical details:",
        *[f"  {key}: {value}" for key, value in wave.data.items()],
    ]
    page.find("#raschii_out").textContent = "\n".join(outputs)
    current_raschii_wave = wave
    page.find("#raschii_plot").innerHTML = wave_profile_to_svg(x, eta)


@pyscript.when("click", "#raschii_plot")
def show_info_when_clicking_plot(mouse_event):
    """
    This is called when the user clicks on the SVG wave profile.

    - Get the coordinates of the clicked point in the wave coordinate system
    - Show information about the clicked point, including particle velocities

    """
    global current_raschii_wave
    if current_raschii_wave is None:
        page.find("#raschii p.info").textContent = "No wave data available."
        return

    # Get the coordinates of the clicked point in the wave coordinate system
    x, z = get_physical_coordinates(mouse_event)

    # Show the information about the clicked point
    info = f"You clicked on x = {x:.3f} m, z = {z:.3f} m (from the bottom)"
    eta = current_raschii_wave.surface_elevation(x=[x], t=0.0)[0]
    if z > eta:
        info += "<br>(Air)"
    else:
        info += "<br>(Water)"
        vel = current_raschii_wave.velocity(x, z, all_points_wet=True)
        info += f"<br>Horizontal particle velocity: {vel[0, 0]:.3f}"
        info += f"<br>Vertical particle velocity:   {vel[0, 1]:.3f}"

    page.find("#raschii p.info").innerHTML = info


class WaveInput:
    def __init__(self):
        self.wave_model_name: str
        self.height: float
        self.depth: float
        self.length: float
        self.order: float

        self.is_ok: bool = True
        self.warning_message: str = ""
        self.error_message: str = ""

    def load_from_page(self):
        """
        Load the wave input parameters from the user input on the web page
        """
        self.wave_model_name = self._get_page_input_value("wave_model", str)
        self.height = self._get_page_input_value("height", float)
        self.depth = self._get_page_input_value("depth", float)
        self.length = self._get_page_input_value("length", float)
        self.order = self._get_page_input_value("order", int)

        if self.height < 0:
            self.is_ok = False
            self.error_message += "Wave height must be non-negative!\n"
        if self.depth <= 0:
            self.is_ok = False
            self.error_message += "Water depth must be non-negative!\n"
        if self.length <= 0:
            self.is_ok = False
            self.error_message += "Wave length must be positive!\n"

    def _get_page_input_value(self, name: str, converter):
        try:
            element = page.find(f"#raschii_{name}")
        except Exception:
            self.is_ok = False
            self.error_message += f"Input element '{name}' not found.\n"
            return None

        try:
            value = element[0].value
        except IndexError:
            self.is_ok = False
            self.error_message += f"Input element '{name}' is empty!\n"
            return None

        try:
            return converter(value)
        except ValueError:
            self.is_ok = False
            self.error_message += f"Invalid value for '{name}': {value!r} is not a number!\n"
            return None


def wave_profile_to_svg(x: list[float], eta: list[float]) -> str:
    """
    Convert the wave profile to SVG format and returns it as string

    We need to handle the fact that SVG y-coordinates are inverted compared to
    the wave coordinate system. The SVG origin is at the top-left corner and is
    positive downwards, while the wave coordinate system has the origin at the
    still water surface and is positive upwards.
    """
    # Plot extents
    xmin = min(x)
    xmax = max(x)
    eta_min = min(eta)
    eta_max = max(eta)
    ymin = max(eta_min - (eta_max - eta_min) * 0.4, 0.0)
    ymax = eta_max + (eta_max - eta_min) * 0.4

    x_plot = [0.0, *x, x[-1]]
    y_plot = [0.0, *eta, 0.0]
    coords = " ".join(f"{x},{y}" for x, y in zip(x_plot, y_plot))

    dx = xmax - xmin
    dy = ymax - ymin
    svg_contents = [
        '<svg version="1.1" preserveAspectRatio="none"',
        f'viewBox="{xmin} {-ymax} {dx} {dy}">',
        '<g transform="scale(1,-1)">',
        f'<path d="M{coords} z" fill="#687dc1"></path>',
        "</g>",
        "</svg>",
    ]
    return "\n".join(svg_contents)


def get_physical_coordinates(mouse_event) -> tuple[float, float]:
    """
    Convert the mouse click coordinates from the SVG element to the physical
    coordinates used when creating the SVG.
    """
    # Coordinates of the click event in client space
    x = mouse_event.clientX
    y = mouse_event.clientY

    # Convert coordinates to the SVG coordinate system
    svg_element = page["#raschii_plot svg"][0]
    cr = svg_element.getBoundingClientRect()
    fx = (x - cr.left) / cr.width
    fy = (y - cr.top) / cr.height

    # Get the SVG view box and calculate the coordinates used when creating the SVG
    # The viewBox is defined as "min-x min-y width height" (in SVG coordinates)
    vb = svg_element.viewBox.baseVal
    x_wave = vb.x + fx * vb.width
    z_wave = vb.y + fy * vb.height
    z_wave *= -1  # Invert the y-coordinate to match the wave coordinate system

    return x_wave, z_wave
