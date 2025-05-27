import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class MobiusStrip:
    def __init__(self, R=5, w=2, n=200):
        """
        Initialize Mobius strip parameters and precompute the meshgrid.

        R: Distance from the center to the strip's centerline (positive float)
        w: Width of the strip (positive float)
        n: Number of sample points along each parameter (integer > 1)

        By Precomputing the meshgrid here avoids redundant computation
        in every method that needs it, making the code cleaner and faster.
        """
        if R <= 0 or w <= 0 or n <= 1:
            raise ValueError("R and w must be positive, and n must be greater than 1.")

        self.R = R
        self.w = w
        self.n = n

        # Precompute parameter ranges
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w / 2, w / 2, n)

        # Meshgrid created once and reused everywhere
        self.U, self.V = np.meshgrid(self.u, self.v)

        # To hold computed 3D coordinates later
        self.X, self.Y, self.Z = None, None, None

    def generate_mesh(self):
        """
        Compute 3D (X, Y, Z) mesh points based on the Möbius strip's parametric equation.

        By separating this into its own method makes it reusable 
        and ensures we don't compute the same values multiple times unnecessarily.
        """
        U, V = self.U, self.V

        self.X = (self.R + V * np.cos(U / 2)) * np.cos(U)
        self.Y = (self.R + V * np.cos(U / 2)) * np.sin(U)
        self.Z = V * np.sin(U / 2)

    def compute_surface_area(self):
        """
        Approximate the surface area using the double integral formula:
        A = ∬ |∂r/∂u × ∂r/∂v| du dv

        Reason: Using discrete finite differences over the meshgrid 
        allows us to numerically approximate the integral with good accuracy.
        """
        du = (2 * np.pi) / (self.n - 1)
        dv = (self.w) / (self.n - 1)

        U, V = self.U, self.V

        # Compute partial derivatives w.r.t u and v
        dx_du = - (self.R + V * np.cos(U / 2)) * np.sin(U) - 0.5 * V * np.sin(U / 2) * np.cos(U)
        dy_du =   (self.R + V * np.cos(U / 2)) * np.cos(U) - 0.5 * V * np.sin(U / 2) * np.sin(U)
        dz_du =  0.5 * V * np.cos(U / 2)

        dx_dv = np.cos(U / 2) * np.cos(U)
        dy_dv = np.cos(U / 2) * np.sin(U)
        dz_dv = np.sin(U / 2)

        # Cross product |∂r/∂u × ∂r/∂v| to get local surface element area
        cross_x = dy_du * dz_dv - dz_du * dy_dv
        cross_y = dz_du * dx_dv - dx_du * dz_dv
        cross_z = dx_du * dy_dv - dy_du * dx_dv

        cross_mag = np.sqrt(cross_x**2 + cross_y**2 + cross_z**2)

        # Double integral approximation via discrete sum
        area = np.sum(cross_mag) * du * dv
        return area

    def compute_edge_length(self):
        """
        Approximate total edge length numerically by summing distances 
        between successive boundary points.

        Since Mobius strips have a continuous single edge, 
        we consider both top and bottom parameter boundaries (v = ±w/2)
        and sum their lengths.

        Numerical integration via discrete point-to-point distances 
        is straightforward and robust for this kind of geometry.
        """
        if self.X is None:
            self.generate_mesh()

        # Extract points along top and bottom edges
        edge_top = np.column_stack((self.X[0, :], self.Y[0, :], self.Z[0, :]))
        edge_bottom = np.column_stack((self.X[-1, :], self.Y[-1, :], self.Z[-1, :]))

        def compute_length(points):
            # Compute Euclidean distances between consecutive points
            diffs = np.diff(points, axis=0)
            segment_lengths = np.linalg.norm(diffs, axis=1)
            total_length = np.sum(segment_lengths)
            # Close the loop by adding distance from last to first point
            total_length += np.linalg.norm(points[-1] - points[0])
            return total_length

        length_top = compute_length(edge_top)
        length_bottom = compute_length(edge_bottom)

        total_length = length_top + length_bottom
        return total_length

    def plot(self, save_path="mobius_strip_plot.png", use_colormap=False):
        """
        Plot the Möbius strip using Matplotlib 3D.

        use_colormap: If True, apply a 'viridis' colormap for better depth perception.
                      Otherwise, use a clean, professional solid color.

        Because a colormap can reveal depth, height, or twist variations more intuitively
        for viewers, while solid colors give a minimal, clean academic-style plot.
        """
        if self.X is None:
            self.generate_mesh()

        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

        if use_colormap:
            # Normalize Z-values for colormap mapping
            norm = (self.Z - self.Z.min()) / (self.Z.max() - self.Z.min())
            ax.plot_surface(self.X, self.Y, self.Z, rstride=5, cstride=5,
                            facecolors=plt.cm.viridis(norm),
                            edgecolor='gray', alpha=0.95)
        else:
            ax.plot_surface(self.X, self.Y, self.Z, rstride=5, cstride=5,
                            color='cornflowerblue', edgecolor='gray', alpha=0.85)

        ax.set_title("3D Möbius Strip")
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.set_zlabel("Z-axis")
        ax.view_init(elev=30, azim=30)

        plt.tight_layout()
        plt.savefig(save_path)
        plt.show()


# Just for example consider we are using unit convention as cm, for our case
# And let's run the code with example parameters for R, w and n.
if __name__ == "__main__":
    strip = MobiusStrip(R=5, w=2, n=200)

    surface_area = strip.compute_surface_area()
    edge_length = strip.compute_edge_length()              #The outputs for the example parameters are

    print(f"Approximate Surface Area: {surface_area:.3f}") #Output: 63.573 cm²
    print(f"Approximate Edge Length: {edge_length:.3f}")   #Output: 67.149 cm

    # Plot using colormap for better 3D depth perception
    strip.plot(use_colormap=True)
