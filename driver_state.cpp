#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    unsigned max = width*height; 
    state.image_color= new pixel[max]; //Sets the array to the resolution of the screen.
    for(size_t i = 0; i < max; i++){
        state.image_color[i] = make_pixel(0, 0, 0); //initializes all colors to black,
    }
    state.image_depth=0;
    std::cout<<"initialize state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    data_geometry *tri;
    unsigned index;
    float *ptr;
    switch(type){
        case render_type::triangle:
            tri = new data_geometry[3];
            index = 0;
            ptr = state.vertex_data;
            for(unsigned i = 0; i <(unsigned) state.num_vertices; i++){
                tri[index].data = ptr;
                index++;
                ptr+=(state.floats_per_vertex-1);
                if(index == 2){
                    index = 0;
                    rasterize_triangle(state,(const data_geometry**) &tri);
                }
            }
            delete []tri;
            break;
        case render_type::indexed:
            
            break;
        case render_type::fan:
            
            break;
        case render_type::strip:
            
            break;
        default:
            break;
    }
    std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    data_vertex v;
    data_geometry g0, g1, g2;
    g0 = (*in)[0]; 
    g1 = (*in)[1];
    g2 = (*in)[2];
    
    v.data = g0.data;
    state.vertex_shader(v, g0, state.uniform_data); 
    v.data = g1.data;
    state.vertex_shader(v, g1, state.uniform_data); 
    v.data = g2.data;
    state.vertex_shader(v, g2, state.uniform_data); 
}

