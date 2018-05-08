def test_readme_example():
    import raschii
    
    fwave = raschii.FentonWave(height=0.25, depth=0.5, length=2.0, order=20)
    print(fwave.coefficients)
    print(fwave.surface_elevation(x=0))
    print(fwave.velocity(x=0, z=0.2))
