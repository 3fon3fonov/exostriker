"""Small SWIFT library example: python exostriker/lib/swift/swift_test.py.

Build first with python exostriker/compile_swift.py.
The automated regression suite is tests/test_swift_library.py.
"""
if __package__:
    from . import integrate
else:
    from api import integrate


def main():
    planets = [
        [0.010956056701692769, 1.286350183286512, 0.12495909431420299,
         90.0, 8.761185986070139, 0.0, 236.93745955638647],
        [0.01207470374287773, 4.267216547236886, 0.19126053191013156,
         90.0, 10.005317383453992, 0.0, 221.47361992077785],
    ]
    result = integrate(1.0, planets, timemax=100.0, timestep=5.0,
                       output_step=1.0, integrator='symba')
    print('Final energy record:', result['energy'][-1])
    for index, orbit in enumerate(result['orbits'], start=1):
        print(f'Planet {index}: {len(orbit)} samples; final orbit:', orbit[-1])


if __name__ == '__main__':
    main()
