import numpy as np
import os
import matplotlib.pyplot as plt
from extra_packages.OutputInterface import OutputInterface
import extra_packages.sperhical_expansion as se


# INPUT STUFF
file_name = 'output_files/CHBrClF1.out'
save_name = 'matched_r.txt'
continue_file = ''

n_plot_points = 50
n_expansion_points = 50
r_plot_list = np.linspace(1, 15, n_plot_points)


# THE PROGRAM

os.system('clear')
print("Let's find some good r values!")
print(f"Loading file: " + file_name.split('/')[-1])

output = OutputInterface(file_name)
Ip = -output.saved_orbitals[output.HOMO][0]  # Should we put this to 1/2?
kappa = np.sqrt(2*Ip)


def radial(r):
    return r**(1/kappa - 1) * np.exp(-kappa*r)


print("Calculating the Laplace coefficients for several r values...")

f_lms = []
for i, r in enumerate(r_plot_list):
    print(f'Evaluating at r={r:.4f} \t Nr. {i+1}/{n_plot_points}')
    f_lms.append(se.spherical_expansion(lambda theta, phi: output.eval_orbital_spherical(r, theta, phi), n_expansion_points))
f_lms = np.array(f_lms)

os.system('clear')
print("Now let's find the good r values!")


max_l = f_lms[0].shape[1]
r_list = []
rl_list = []
l = 0
m = 0

if continue_file:  # Load the values we need to continue on
    print('You have choosen to contiune on ' + continue_file)
    r_list = list(np.loadtxt(continue_file))
    if np.sqrt(len(r_list)) % 1 != 0:
        raise ValueError('The length of the file to continue on does not match a l value!')
    l = int(np.sqrt(len(r_list)))
    m = -l

while True:
    if l == max_l:
        print("You have reached the max l in the Laplace expansion! Stopping and saving...")
        break

    print(f"Currently for : l = {l} and m = {m}")
    sign = 0 if m >= 0 else 1

    # Make the plot
    plt.figure(facecolor='white')
    plt.plot(r_plot_list, np.abs(f_lms[:, sign, l, abs(m)]) * 100, label='Laplace expan. coeffs.')
    plt.plot(r_plot_list, np.abs(f_lms[:, sign, l, abs(m)]) / radial(r_plot_list), label='Asymptotic coeff.')
    plt.xlabel('$r$ (a.u.)')
    plt.ylabel('Absolute amplitude')
    plt.legend()
    plt.show()

    user_input = input('Input : ')

    if user_input == 'exit':
        break
    else:
        try:
            val = float(user_input)
        except ValueError:
            print('Not a number or unknown command!')
            continue  # Try again

    rl_list.append(val)

    # Continue to next l if m is max value
    if m == l:
        l += 1
        m = -l
        r_list = r_list + rl_list  # Add the values for the current l to the end
        rl_list = []
    else:
        m = m + 1

    os.system('clear')

# Now save the list!
np.savetxt(save_name, r_list)


