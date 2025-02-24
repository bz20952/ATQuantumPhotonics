def rect_beam(b, h):

    A = b*h
    I_xx = b*h**3/12

    return I_xx, A


def t_beam(H, b, h, B):

    """See equation at https://amesweb.info/section/second-moment-of-area-calculator.aspx."""

    A = B*h + b*H
    y_c = ((H+h/2)*h*B + H**2*b/2)/A
    I_xx = b*H*(y_c-H/2)**2 + b*H**3/12 + h*B*(H+h/2-y_c)**2 + h**3*B/12

    return I_xx, A


def i_beam(H, b, h, B):

    """See equation at https://amesweb.info/section/second-moment-of-area-calculator.aspx."""

    A = 2*h*B + b*H
    I_xx = H**3*b/12 + 2*(h**3*B/12 + h*B*(H+h)**2/4)
    
    return I_xx, A