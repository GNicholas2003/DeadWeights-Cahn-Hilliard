        import numpy as np
        import matplotlib.pyplot as plt
        
        # ----------------------------------------------------
        # Parameters
        # ----------------------------------------------------
        N = 128              # grid size
        dx = 1.0
        dt = 0.01
        
        M = 1.0              # mobility
        kappa = 1.0
        
        steps = 5000
        plot_every = 500
        
        # ----------------------------------------------------
        # Initial condition
        # ----------------------------------------------------
        c = 0.05*np.random.randn(N, N)
        
        # ----------------------------------------------------
        # Periodic Laplacian
        # ----------------------------------------------------
        def laplacian(u):
            return (
                np.roll(u,1,0)+np.roll(u,-1,0)+
                np.roll(u,1,1)+np.roll(u,-1,1)
                -4*u
            )/dx**2
        
        # ----------------------------------------------------
        # Time integration
        # ----------------------------------------------------
        plt.figure(figsize=(6,6))
        
        for step in range(steps):
        
            lap_c = laplacian(c)
        
            # chemical potential
            mu = c**3 - c - kappa*lap_c
        
            # Cahn-Hilliard equation
            c += dt * M * laplacian(mu)
        
            if step % plot_every == 0:
                plt.clf()
                plt.imshow(c,
                           cmap="RdBu",
                           origin="lower",
                           vmin=-1,
                           vmax=1)
                plt.title(f"Step {step}")
                plt.colorbar()
                plt.pause(0.01)
        
        plt.show()
