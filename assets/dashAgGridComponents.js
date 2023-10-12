var dagcomponentfuncs = window.dashAgGridComponentFunctions = window.dashAgGridComponentFunctions || {};
var dagfuncs = window.dashAgGridFunctions = window.dashAgGridFunctions || {};

function isHeaderRow(params) {
    return params.data.section === 'big-title';
}

dagfuncs.titleColSpan = function (params, ncols) {
    if (isHeaderRow(params)) {
        return ncols;
    } else {
        return 1;
    }
}

//Button
dagcomponentfuncs.Button = function (props) {
    const {setData, data} = props;

    function onClick() {
        setData();
    }
    return React.createElement(
        'button',
        {
            onClick: onClick,
            className: props.className,
        },
        props.value
    );
};

//Dropdown

dagfuncs.dynamicOptions = function(params, values) {
    return {
            values: values
        };
    }
dagfuncs.DCC_Dropdown = class {

    // gets called once before the renderer is used
  init(params) {
    // create the cell
    this.params = params;
    this.ref = React.createRef();

    // function for when Dash is trying to send props back to the component / server
    var setProps = (props) => {
        if (props.value) {
            // updates the value of the editor
            this.value = props.value;

            // re-enables keyboard event
            delete params.colDef.suppressKeyboardEvent

            // tells the grid to stop editing the cell
            params.api.stopEditing();

            // sets focus back to the grid's previously active cell
            this.prevFocus.focus();
        }
    }
    this.eInput = document.createElement('div')

    // renders component into the editor element
    ReactDOM.render(React.createElement(window.dash_core_components.Dropdown, {
        options: params.values, value: params.value, clearable: true, ref: this.ref, setProps, style: {width: params.column.actualWidth},
    }), this.eInput)

    // allows focus event
    this.eInput.tabIndex = "0"

    // sets editor value to the value from the cell
    this.value = params.value;
  }

  // gets called once when grid ready to insert the element
  getGui() {
    return this.eInput;
  }

  focusChild() {
    // mousedown event
    const clickEvent = new MouseEvent('mousedown', {
        view: window,
        bubbles: true
    });

    // needed to delay and allow the component to render
    setTimeout(() => {
        var inp = this.eInput.getElementsByClassName('Select-control')[0]

        // disables keyboard event
        this.params.colDef.suppressKeyboardEvent = (params) => {
               const gridShouldDoNothing = params.editing
               return gridShouldDoNothing;
           }
        // shows dropdown options
        inp.dispatchEvent(clickEvent)
    }, 100)
  }

  // focus and select can be done after the gui is attached
  afterGuiAttached() {
    // stores the active cell
    this.prevFocus = document.activeElement

    // adds event listener to trigger event to go into dash component
    this.eInput.addEventListener('focus', this.focusChild())

    // triggers focus event
    this.eInput.focus();
  }

  // returns the new value after editing
  getValue() {
    return this.value;
  }

  // any cleanup we need to be done here
  destroy() {
    // sets focus back to the grid's previously active cell
    this.prevFocus.focus();
  }
}