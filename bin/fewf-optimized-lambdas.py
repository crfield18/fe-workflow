import numpy as np
from scipy.interpolate import interp1d

def s2(x):
    """The second order smoothstep function. 
    
    Parameters
    ----------
    x : float or array_like
        The input values.
    
    Returns
    -------
    float or array_like
        The output values.
    
    """
    return 6*x**5 -15*x**4 +10*x**3

def inverse_s2_interpolation(y_values, x_min=0, x_max=1, kind='cubic'):
    """ Calculate the inverse of the second order smoothstep function using interpolation.
    
    Parameters
    ----------
    y_values : array_like
        The output values.
    x_min : float, optional
        The minimum input value. Default is 0.
    x_max : float, optional
        The maximum input value. Default is 1.
    kind : str, optional
        The kind of interpolation. Default is 'cubic'.

    Returns
    -------
    array_like
        The input values.

    """
    x_interp = np.linspace(x_min, x_max, 1000)  # More points = better accuracy
    y_interp = s2(x_interp)
    f_inverse = interp1d(y_interp, x_interp, kind=kind, bounds_error=False, fill_value=np.nan)
    # Calculate inverse values using interpolation
    return f_inverse(y_values)

def Predict_Lambda(nlambda):
    """ Predict lambda values using the second order smoothstep function.
    
    Parameters
    ----------
    nlambda : int
        The number of lambda values to predict.
    
    Returns
    -------
    array_like
        The predicted lambda values.

    """
    xnew = np.linspace(0,1, nlambda)
    lams = inverse_s2_interpolation(xnew)
    return lams

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Predict lambda values')
    parser.add_argument('--nlambda', type=int, help='Number of lambda values')
    parser.add_argument('--output', type=str, help='Output file')
    args = parser.parse_args()

    # Prints output to the screen
    if args.nlambda:
        lambdas = Predict_Lambda(args.nlambda)
        print(f"lams=({' '.join([f'{x:.8f}' for x in lambdas])})")

        # Writes output to a file
        if args.output:
            with open(args.output, 'w') as f:
                for lam in lambdas:
                    f.write(f"{lam:.8f}\n")
